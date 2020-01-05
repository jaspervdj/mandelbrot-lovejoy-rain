{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE DeriveFunctor       #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main
    ( main
    ) where

import qualified Codec.Picture.Gif                 as JP
import qualified Codec.Picture.Png                 as JP
import qualified Codec.Picture.Types               as JP
import qualified Control.Concurrent.MVar           as MVar
import           Control.Monad                     (foldM, forM, forM_,
                                                    replicateM_, when)
import           Control.Monad.Primitive           (PrimMonad (..))
import           Control.Monad.ST                  (runST)
import qualified Control.Scheduler                 as Scheduler
import           Data.Bool                         (bool)
import qualified Data.IORef                        as IORef
import           Data.List                         (foldl')
import qualified Data.Massiv.Array                 as M
import qualified Data.Massiv.Array.Manifest.Vector as MV
import qualified Data.Massiv.Array.Mutable         as MM
import qualified Data.Massiv.Array.Unsafe          as M
import           Data.Maybe                        (fromMaybe, listToMaybe)
import           Data.Proxy                        (Proxy (..))
import qualified Data.Vector                       as V
import           Data.Word                         (Word32)
import           Data.Word                         (Word8)
import qualified System.IO                         as IO
import qualified System.Random.MWC                 as MWC

--------------------------------------------------------------------------------
-- Utility types.

newtype Radius   = Radius   {unRadius   :: Float}
newtype Distance = Distance Float


--------------------------------------------------------------------------------
-- Utilities for dimension-generic algorithms.

unit :: M.Index ix => ix
unit = M.pureIndex 1

(.-.), (.+.) :: M.Index ix => ix -> ix -> ix
(.-.) = M.liftIndex2 (-)
(.+.) = M.liftIndex2 (+)

intersection :: M.Index ix => (ix, ix) -> (ix, ix) -> (ix, ix)
intersection (l1, r1) (l2, r2) =
    (M.liftIndex2 max l1 l2, M.liftIndex2 min r1 r2)

offset :: M.Index ix => ix -> HPulse ix -> HPulse ix
offset o hp = hp {hpCenter = hpCenter hp .-. o}

distance :: M.Index ix => ix -> ix -> Distance
distance i j =
    Distance . sqrt . fromIntegral . M.foldlIndex (+) 0 $
    M.liftIndex2 (\p s -> (p - s) * (p - s)) i j


--------------------------------------------------------------------------------
-- Generating pulses

data HPulse ix = HPulse
    { hpCenter :: !ix
    , hpAmp    :: !Float
    , hpRadius :: !Radius
    }

arbitraryIndex
    :: forall m ix. (PrimMonad m, M.MonadThrow m, M.Index ix)
    => MWC.Gen (PrimState m) -> M.Sz ix -> m ix
arbitraryIndex gen (M.Sz choice) =
    foldM
        (\acc dim -> do
            up <- M.getDimM choice (M.Dim dim)
            x  <- MWC.uniformR (0, up) gen
            M.setDimM acc (M.Dim dim) x)
        choice
        [1 .. M.unDim (M.dimensions (Proxy :: Proxy ix))]

arbitraryPulse
    :: forall m ix. (PrimMonad m, M.MonadThrow m, M.Index ix)
    => MWC.Gen (PrimState m) -> Float -> M.Sz ix -> m (HPulse ix)
arbitraryPulse gen alpha area = do
    center <- arbitraryIndex gen area
    rhoInv <- MWC.uniformR (0.0, 1.0) gen
    let rho   =  if rhoInv == 0.0 then 1.0 else 1.0 / (1.0 - rhoInv)
        radius = rho / 2.0
        amp    = rho ** (1.0 / alpha)
    ampSign <- bool (-1.0) 1.0 <$> MWC.uniform gen
    return $ HPulse center (ampSign * amp) (Radius radius)


--------------------------------------------------------------------------------
-- Drawing pulses

data PulseShape
    = Rectangular
    | Smooth Float    -- s
    | Annuli Float Float  -- lambda, s

pulseAt :: PulseShape -> Radius -> Distance -> Float
pulseAt Rectangular       (Radius r) (Distance u) = if u <= r then 1.0 else 0.0
pulseAt (Smooth s)        (Radius r) (Distance u) = exp (-(u / r) ** (2 * s))
pulseAt (Annuli lambda s) (Radius r) (Distance u) =
    let lambda' = sqrt (lambda * lambda - 1.0)
        delta   = (lambda + lambda') / 2.0
        sigma   = (lambda - lambda') / 2.0 in
    exp (-(((u*u) / (r*r) - (delta*delta)) / (sigma*sigma)) ** (2 * s))

bounds :: M.Index ix => PulseShape -> HPulse ix -> (ix, ix)
bounds shape hp =
    let bound = fromMaybe 0.0 . listToMaybe .
            dropWhile ((>= 0.1e-6) . pulseAt shape (hpRadius hp) . Distance) .
            iterate succ . unRadius $ hpRadius hp

        radIdx = M.pureIndex $ ceiling bound
        start  = hpCenter hp .-. radIdx
        end    = hpCenter hp .+. radIdx .+. unit in
    (start, end)

addPulse
    :: forall ix r m. (MM.Mutable r ix Float, M.PrimMonad m, M.MonadThrow m)
    => MM.MArray (M.PrimState m) r ix Float
    -> PulseShape
    -> HPulse ix
    -> m ()
addPulse marr shape pulse@HPulse {..} = M.iterM_ i0 end unit (<) $ \i ->
    MM.modifyM
        marr
        (\x -> pure $ x + hpAmp * pulseAt shape hpRadius (distance i hpCenter))
        i
  where
    pulseBounds = bounds shape pulse
    marrBounds  = (M.pureIndex 0, M.unSz (M.msize marr))
    (i0, end)   = intersection pulseBounds marrBounds
{-# INLINE addPulse #-}


--------------------------------------------------------------------------------
-- Generating the fractal.

data FractalParams ix = FractalParams
    { fpWorld      :: M.Sz ix
    , fpCropOffset :: ix
    , fpCropSize   :: M.Sz ix
    , fpNumPulses  :: Int
    , fpAlpha      :: Float
    , fpPulseShape :: PulseShape
    , fpRhoIn      :: Float
    , fpRhoOut     :: Float
    , fpTreshold   :: Float
    }

makeFractal
    :: (M.Index ix, PrimMonad m, M.MonadThrow m)
    => FractalParams ix
    -> MWC.Gen (PrimState m)
    -> m (M.Array M.P ix Float)
makeFractal FractalParams {..} gen = do
    marr <- MM.new fpCropSize
    replicateM_ fpNumPulses $ do
        pulse <- offset fpCropOffset <$> arbitraryPulse gen fpAlpha fpWorld
        when (relevant pulse) $ addPulse marr fpPulseShape pulse
    M.unsafeFreeze M.Seq marr
  where
    relevant = (\r -> r >= fpRhoIn && r <= fpRhoOut) . unRadius . hpRadius

--------------------------------------------------------------------------------
-- | Post-processing utilities.

normalize
    :: (Functor (M.Array r ix), M.Source r ix Float)
    => M.Array r ix Float -> M.Array r ix Float
normalize arr =
    let (mini, maxi) = (M.minimum' arr, M.maximum' arr) in
    (\x -> (x - mini) / (maxi - mini)) <$> arr

treshold
    :: forall r ix. (Functor (M.Array r ix), M.Source r ix Float)
    => Float -> M.Array r ix Float -> M.Array r ix Float
treshold relativeTreshold arr =
    (\x ->
        if x < absoluteTreshold
            then 0.0
            else (x - absoluteTreshold) / (1.0 - absoluteTreshold)) <$> arr
  where
    numBuckets = 100

    histogram :: M.Array M.P Int Int
    histogram = runST $ do
        marr <- MM.makeMArrayS (M.Sz numBuckets + 1) (\_ -> pure 0)
        M.forM_ arr $ \x ->
            let b = floor (x * fromIntegral numBuckets) in
            MM.modifyM marr (pure . succ) b
        MM.freezeS marr

    absoluteTreshold =
        let target = floor $
                relativeTreshold * fromIntegral (M.totalElem (M.size arr))
            go acc i
                | i >= numBuckets = 1.0
                | otherwise       =
                    let acc' = acc + fromMaybe 0 (M.index histogram i) in
                    if acc' >= target
                        then fromIntegral i / fromIntegral numBuckets
                        else go acc' (i + 1)  in
        go 0 0

floatToWord8 :: Float -> Word8
floatToWord8 = round . (* 255)

pixelToWord8 :: JP.PixelRGBF -> JP.PixelRGB8
pixelToWord8 (JP.PixelRGBF r g b) =
    JP.PixelRGB8 (floatToWord8 r) (floatToWord8 g) (floatToWord8 b)

palette :: JP.Palette
palette = JP.generateImage
    (\x _ -> pixelToWord8 $ interpolate (fromIntegral x / 255.0) blue white)
    256
    1

dark, blue, white :: JP.PixelRGBF
dark  = JP.PixelRGBF 0.0 0.1 0.3
blue  = JP.PixelRGBF 0.0 0.4 1.0
white = JP.PixelRGBF 1.0 1.0 1.0

class Interpolate a where
    interpolate :: Float -> a -> a -> a

instance Interpolate Float where
    interpolate t x y = (1.0 - t) * x + t * y

instance Interpolate JP.PixelRGBF where
    interpolate t (JP.PixelRGBF r0 g0 b0) (JP.PixelRGBF r1 g1 b1) = JP.PixelRGBF
        (interpolate t r0 r1) (interpolate t g0 g1) (interpolate t b0 b1)

interpolateFrames
    :: M.Index ix
    => Int  -- Number of interpolation frames to insert
    -> [M.Array M.D ix Float]
    -> [M.Array M.D ix Float]
interpolateFrames _   []           = []
interpolateFrames _   [x]          = [x]
interpolateFrames num (x : y : ys) =
    [x] ++
    [ M.zipWith (interpolate (fromIntegral n * spacing)) x y
    | n <- [1 .. num]
    ] ++
    interpolateFrames num (y : ys)
  where
    spacing = 1.0 / fromIntegral (num + 1)


--------------------------------------------------------------------------------
-- | Parallellizing the drawing across cores

mapReduce
    :: (String -> IO ())  -- ^ Logger
    -> [a]                -- ^ Items to process
    -> (a -> b)           -- ^ Map
    -> (b -> b -> b)      -- ^ Reduce
    -> IO b
mapReduce logger as f g = do
    bDone <- IORef.newIORef (0 :: Int)
    bs    <-
        Scheduler.withScheduler Scheduler.Par $ \scheduler ->
            forM_ as $ \a ->
            Scheduler.scheduleWork scheduler $ do
            let b = f a
            done <- b `seq` IORef.atomicModifyIORef' bDone $ \n -> (succ n, succ n)
            logger $ "Finished job " ++ show done ++ "/" ++ show num
            pure b
    case bs of
        []      -> fail "mapReduce: empty list"
        bh : bt -> return $ foldl' g bh bt
  where
    num = length as

divideWork :: Int -> Int -> [Int]
divideWork total workers =
    let (items, remainder) = total `divMod` workers in
    zipWith (+) (replicate remainder 1 ++ repeat 0) (replicate workers items)

makeClouds
    :: M.Index ix
    => (String -> IO ())   -- ^ Logging
    -> FractalParams ix    -- ^ Parameters
    -> IO (M.Array M.P ix Float)
makeClouds logger fp@FractalParams {..} = do
    jobs <- forM (divideWork fpNumPulses 32) $ \items -> do
        seed <- MWC.withSystemRandom $ MWC.asGenIO $
            V.replicateM 32 . MWC.uniform
        return (items :: Int, seed :: V.Vector Word32)

    summed <- mapReduce logger jobs
        (\(items, seed) -> runST $ do
            gen <- MWC.initialize seed
            makeFractal fp {fpNumPulses = items} gen)
        (\x y -> M.computeProxy (Proxy :: Proxy M.P) $ M.zipWith (+) x y)

    return $ M.computeProxy (Proxy :: Proxy M.P) $
        fmap scurve $ treshold fpTreshold $
        normalize $ M.delay summed
  where
    scurve x = (1 + cos ((x - 1) * pi)) / 2

makeGradientBackground
    :: PrimMonad m
    => MWC.Gen (PrimState m) -> Int -> Int -> m (JP.Image JP.PixelRGBF)
makeGradientBackground gen width height = do
    phi <- MWC.uniformR (0, 2 * pi) gen
    let measure = fromIntegral $ max width height
        rho     = 3 * measure
        cx      = fromIntegral width  * 0.5 + rho * cos phi
        cy      = fromIntegral height * 0.5 + rho * sin phi
        render  = \x y ->
            let dx   = fromIntegral x - cx
                dy   = fromIntegral y - cy
                dist = sqrt $ dx * dx + dy * dy
                t    = (dist - (rho - measure)) / (2 * measure) in
            interpolate t blue dark
    pure $ JP.generateImage render width height

main2d :: IO ()
main2d = do
    let height   =  900
        width    = 1600
        small    = M.Ix2 height width
        world    = M.Sz (small .+. small .+. small)
        v        = 3

        fp = FractalParams
            { fpWorld      = world
            , fpCropOffset = small
            , fpCropSize   = M.Sz small
            , fpNumPulses  = M.totalElem world * v
            , fpAlpha      = 5 / 3
            , fpPulseShape = Rectangular -- Annuli 1.2 8
            , fpRhoIn      = 0.0
            , fpRhoOut     = fromIntegral (max height width) * 1.5
            , fpTreshold   = 0.6
            }

    -- Logging
    lock <- MVar.newMVar ()
    let logger = MVar.modifyMVar_ lock . const . IO.hPutStrLn IO.stderr

    gradientBackground <- MWC.withSystemRandom $ \gen ->
        (makeGradientBackground gen width height :: IO (JP.Image JP.PixelRGBF))

    -- Blend the background and the foreground.
    clouds <- fmap
        (MV.toVector .  M.computeProxy (Proxy :: Proxy M.P) . M.delay)
        (makeClouds logger fp)
    let image :: JP.Image JP.PixelRGBF
        image = JP.generateImage
            (\x y ->
                let bg    = JP.pixelAt gradientBackground x y
                    cloud = clouds V.! (x + y * width) in
                interpolate cloud bg white)
            width
            height

    JP.writePng "2d.png" $ JP.pixelMap pixelToWord8 image

main3d :: IO ()
main3d = do
    let numkeyframes =  32
        nonkeyframes =   2
        delay        =   4

        height   = 200
        width    = 300
        small    = M.Ix3 numkeyframes height width
        world    = M.Sz (1 .+. small .+. small)
        v        = 1

        fp = FractalParams
            { fpWorld      = world
            , fpCropOffset = small
            , fpCropSize   = M.Sz small
            , fpNumPulses  = M.totalElem world * v `div` 20
            , fpAlpha      = 5 / 3
            , fpPulseShape = Annuli 1.2 8
            , fpRhoIn      = 0.0
            , fpRhoOut     = fromIntegral (max height width) * 1.5
            , fpTreshold   = 0.6
            }

    -- Logging
    lock <- MVar.newMVar ()
    let logger = MVar.modifyMVar_ lock . const . IO.hPutStrLn IO.stderr

    cuby <- makeClouds logger fp

    logger "Generating keyframes..."
    let keyframes :: [M.Array M.M M.Ix2 Float]
        keyframes = [cuby M.!> d | d <- [0 .. numkeyframes - 1]]

        frames :: [JP.Image Word8]
        frames =
            [ JP.Image width height $ MV.toVector $
                M.computeProxy (Proxy :: Proxy M.P) $
                fmap floatToWord8 frame
            | frame <- interpolateFrames nonkeyframes (map M.delay keyframes)
            ]

    either fail id $ JP.writeComplexGifImage "3d.gif" JP.GifEncode
        { JP.geWidth      = width
        , JP.geHeight     = height
        , JP.gePalette    = Just palette
        , JP.geBackground = Nothing
        , JP.geLooping    = JP.LoopingForever
        , JP.geFrames     = do
            frame <- frames
            return JP.GifFrame
                { JP.gfXOffset     = 0
                , JP.gfYOffset     = 0
                , JP.gfPalette     = Nothing
                , JP.gfTransparent = Nothing
                , JP.gfDelay       = delay
                , JP.gfDisposal    = JP.DisposalDoNot
                , JP.gfPixels      = frame
                }
        }

main :: IO ()
main = main2d >> main3d
