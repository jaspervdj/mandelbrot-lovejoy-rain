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
import qualified Codec.Picture.Types               as JP
import           Control.Applicative
import           Control.Concurrent                (forkIO)
import qualified Control.Concurrent.Chan           as Chan
import qualified Control.Concurrent.MVar           as MVar
import           Control.Monad                     (foldM, forM, forM_,
                                                    replicateM_, when)
import           Control.Monad.Primitive           (PrimMonad (..))
import           Control.Monad.ST                  (runST)
import qualified Control.Scheduler                 as Scheduler
import           Data.Bool                         (bool)
import qualified Data.IORef                        as IORef
import           Data.List                         (foldl')
import           Data.Massiv.Array                 as M
import qualified Data.Massiv.Array.Mutable         as MM
import qualified Data.Massiv.Array.Unsafe          as M
import qualified Data.Massiv.Array.IO              as MIO
import           Data.Maybe                        (catMaybes, fromMaybe,
                                                    isJust, listToMaybe)
import qualified Data.Vector                       as V
import           Data.Word                         (Word32, Word8)
import           Graphics.ColorSpace
import qualified System.IO                         as IO
import qualified System.Random.MWC                 as MWC

--------------------------------------------------------------------------------
-- Utility types.

newtype Radius   = Radius   {unRadius   :: Float}
newtype Distance = Distance Float


--------------------------------------------------------------------------------
-- Utilities for dimension-generic algorithms.

intersection :: M.Index ix => (ix, ix) -> (ix, ix) -> (ix, ix)
intersection (l1, r1) (l2, r2) =
    (M.liftIndex2 max l1 l2, M.liftIndex2 min r1 r2)

offset :: (Num ix, M.Index ix) => ix -> HPulse ix -> HPulse ix
offset o hp = hp {hpCenter = hpCenter hp - o}

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
arbitraryIndex gen sz@(M.Sz choice) =
    foldM
        (\acc dim -> do
            up <- M.getDimM choice dim
            x  <- MWC.uniformR (0, up) gen
            M.setDimM acc dim x)
        choice
        [1 .. M.dimensions sz]

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

bounds :: (Num ix, M.Index ix) => PulseShape -> HPulse ix -> (ix, ix)
bounds shape hp =
    let bound = fromMaybe 0.0 . listToMaybe .
            dropWhile ((>= 0.1e-6) . pulseAt shape (hpRadius hp) . Distance) .
            iterate succ . unRadius $ hpRadius hp

        radIdx = M.pureIndex $ ceiling bound
        start  = hpCenter hp - radIdx
        end    = hpCenter hp + radIdx + M.oneIndex in
    (start, end)

addPulse
    :: forall ix r m. (Num ix, MM.Mutable r ix Float, M.PrimMonad m, M.MonadThrow m)
    => MM.MArray (M.PrimState m) r ix Float
    -> PulseShape
    -> HPulse ix
    -> m ()
addPulse marr shape pulse@HPulse {..} = M.iterM_ i0 end M.oneIndex (<) $ \i ->
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
    { fpWorld      :: !(M.Sz ix)
    , fpCropOffset :: !ix
    , fpCropSize   :: !(M.Sz ix)
    , fpNumPulses  :: !Int
    , fpAlpha      :: !Float
    , fpPulseShape :: !PulseShape
    , fpRhoIn      :: !Float
    , fpRhoOut     :: !Float
    , fpThreshold  :: !Float
    }


-- writeFractal ::
--        (Num ix, M.Index ix, PrimMonad m, M.MonadThrow m)
--     => FractalParams ix
--     -> MWC.Gen (PrimState m)
--     -> M.MArray (PrimState m) M.P ix Float
--     -> m ()
-- writeFractal FractalParams {..} gen marr = do
--     replicateM_ fpNumPulses $ do
--         pulse <- offset fpCropOffset <$> arbitraryPulse gen fpAlpha fpWorld
--         when (relevant pulse) $ addPulse marr fpPulseShape pulse
--   where
--     relevant = (\r -> r >= fpRhoIn && r <= fpRhoOut) . unRadius . hpRadius

makeFractal
    :: (Num ix, M.Index ix, PrimMonad m, M.MonadThrow m)
    => FractalParams ix
    -> MWC.Gen (PrimState m)
    -> m (M.Array M.P ix Float)
makeFractal FractalParams {..} gen =
    MM.createArrayS_ M.Seq fpCropSize $ \ marr -> replicateM_ fpNumPulses $ do
        pulse <- offset fpCropOffset <$> arbitraryPulse gen fpAlpha fpWorld
        when (relevant pulse) $ addPulse marr fpPulseShape pulse
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

threshold
    :: forall r ix. (Functor (M.Array r ix), M.Source r ix Float)
    => Float -> M.Array r ix Float -> M.Array r ix Float
threshold relativeThreshold arr =
    (\x ->
        if x < absoluteThreshold
            then 0.0
            else (x - absoluteThreshold) / (1.0 - absoluteThreshold)) <$> arr
  where
    numBuckets = 100

    histogram :: M.Array M.P Int Int
    histogram = runST $ do
        marr <- MM.makeMArrayS (M.Sz numBuckets + 1) (\_ -> pure 0)
        M.forM_ arr $ \x ->
            let b = floor (x * fromIntegral numBuckets) in
            MM.modifyM marr (pure . succ) b
        MM.freezeS marr

    absoluteThreshold =
        let target = floor $
                relativeThreshold * fromIntegral (M.totalElem (M.size arr))
            go acc i
                | i >= numBuckets = 1.0
                | otherwise       =
                    let acc' = acc + fromMaybe 0 (M.index histogram i) in
                    if acc' >= target
                        then fromIntegral i / fromIntegral numBuckets
                        else go acc' (i + 1)  in
        go 0 0

palette :: JP.Palette
palette =
    MIO.toJPImageRGB8 $
    M.resize' (Sz2 1 256) $
    M.makeArrayR D Seq (Sz1 256) $ \x ->
      toWord8 $ interpolate (fromIntegral x / 255.0) blue white

dark, blue, white :: Pixel RGB Float
dark  = PixelRGB 0.0 0.1 0.3
blue  = PixelRGB 0.0 0.4 1.0
white = PixelRGB 1.0 1.0 1.0

class Interpolate a where
    interpolate :: Float -> a -> a -> a

instance Interpolate Float where
    interpolate t x y = (1.0 - t) * x + t * y

instance Interpolate (Pixel RGB Float) where
    interpolate t = liftA2 (interpolate t)

interpolateFrames
    :: Int  -- Number of interpolation frames to insert
    -> M.Array M.B M.Ix1 (M.Array M.M Ix2 Float)
    -> M.Array M.B M.Ix1 (MIO.Image M.S Y Word8)
interpolateFrames num arr =
    maybe (M.compute (M.map toArrayY arr)) (M.compute . setComp Par) $ do
      (_, arr') <- M.unconsM arr
      (_, l) <- M.unsnocM arr
      frames <- concatM 1 $ M.zipWith interpolateSingle arr arr'
      pure $ M.snoc frames $ toArrayY l
  where
    spacing = 1.0 / fromIntegral (num + 1)
    toArrayY :: Source r Ix2 Float => M.Array r Ix2 Float -> MIO.Image M.S Y Word8
    toArrayY = M.compute . M.map (PixelY . eToWord8)
    interpolateSingle ::
      M.Array M.M Ix2 Float -> M.Array M.M Ix2 Float -> M.Array M.D Ix1 (MIO.Image M.S Y Word8)
    interpolateSingle x y =
      (\n -> toArrayY (M.zipWith (interpolate (fromIntegral n * spacing :: Float)) x y)) <$>
      (1 ... num)


--------------------------------------------------------------------------------
-- | Parallellizing the drawing across cores

divideWork :: Int -> Int -> [Int]
divideWork total workers =
    let (items, remainder) = total `divMod` workers in
    Prelude.zipWith (+) (Prelude.replicate remainder 1 ++ repeat 0) (Prelude.replicate workers items)

weights :: Int -> Int -> M.Array M.P Int Int
weights fpNumPulses nWorkers = M.fromList Par $ divideWork fpNumPulses nWorkers


-- extractWorkerStates = _workerStatesArray

-- makeClouds
--     :: (Num ix, M.Index ix)
--     => (String -> IO ())   -- ^ Logging
--     -> FractalParams ix    -- ^ Parameters
--     -> IO (M.Array M.P ix Float)
-- makeClouds logger fp@FractalParams {..} = do
--     wss <- initWorkerStates Par $ \(ThreadId i) -> do
--         gen <- MWC.createSystemRandom
--         marr <- MM.new fpCropSize
--         pure (i, gen, marr)
--     void $ M.forWS wss (weights fpNumPulses 32) $ \items (_, gen, marr) ->
--         writeFractal fp {fpNumPulses = items} gen marr
--     withSchedulerWS wss $ \scheduler -> unwrapSchedulerWS scheduler
--     jobs <- Control.Monad.forM (divideWork fpNumPulses 32) $ \items -> do
--         seed <- MWC.withSystemRandom $ MWC.asGenIO $
--             V.replicateM 32 . MWC.uniform
--         return (items :: Int, seed :: V.Vector Word32)

--     summed <- mapReduce logger jobs
--         (\(items, seed) -> runST $ do
--             gen <- MWC.initialize seed
--             makeFractal fp {fpNumPulses = items} gen)
--         (\x y -> M.computeAs M.P $ M.zipWith (+) x y)

--     return $ M.computeAs M.P $
--         fmap scurve $ threshold fpThreshold $
--         normalize $ M.delay summed
--   where
--     scurve x = (1 + cos ((x - 1) * pi)) / 2


makeClouds
    :: (Num ix, M.Index ix)
    => (String -> IO ())   -- ^ Logging
    -> FractalParams ix    -- ^ Parameters
    -> IO (M.Array M.P ix Float)
makeClouds logger fp@FractalParams {..} = do
    bDone <- IORef.newIORef (0 :: Int)
    gens <- initWorkerStates Par (const MWC.createSystemRandom)
    let jobs = weights fpNumPulses 32
        num = M.unSz (M.size jobs)
    parClouds :: M.Array M.B Ix1 (M.Array M.P ix Float) <-
        M.forWS gens jobs $ \items gen -> do
            arr <- makeFractal fp {fpNumPulses = items} gen
            done <- IORef.atomicModifyIORef' bDone $ \n -> (succ n, succ n)
            logger $ "Finished job " ++ show done ++ "/" ++ show num
            pure arr
    return $ M.computeAs M.P $
        fmap scurve $ threshold fpThreshold $
        normalize $ M.delay $ foldl1 (\x -> compute . M.zipWith (+) x) parClouds
  where
    scurve x = (1 + cos ((x - 1) * pi)) / 2

makeGradientBackground
    :: PrimMonad m
    => MWC.Gen (PrimState m) -> M.Sz2 -> m (MIO.Image M.D RGB Float)
makeGradientBackground gen sz@(M.Sz2 height width) = do
    phi <- MWC.uniformR (0, 2 * pi) gen
    let measure = fromIntegral $ max width height
        rho     = 3 * measure
        cx      = fromIntegral width  * 0.5 + rho * cos phi
        cy      = fromIntegral height * 0.5 + rho * sin phi
    pure $ M.makeArray Par sz $ \ (y :. x) ->
            let dx   = fromIntegral x - cx
                dy   = fromIntegral y - cy
                dist = sqrt $ dx * dx + dy * dy
                t    = (dist - (rho - measure)) / (2 * measure) in
            interpolate t blue dark


main2d :: IO ()
main2d = do
    let height   =  900
        width    = 1600
        small    = M.Sz2 height width
        world    = 3 * small
        v        = 3

        fp = FractalParams
            { fpWorld      = world
            , fpCropOffset = M.unSz small
            , fpCropSize   = small
            , fpNumPulses  = M.totalElem world * v
            , fpAlpha      = 5 / 3
            , fpPulseShape = Rectangular -- Annuli 1.2 8
            , fpRhoIn      = 0.0
            , fpRhoOut     = fromIntegral (max height width) * 1.5
            , fpThreshold   = 0.6
            }

    -- Logging
    lock <- MVar.newMVar ()
    let logger = MVar.modifyMVar_ lock . const . IO.hPutStrLn IO.stderr

    gradientBackground <- MWC.withSystemRandom . MWC.asGenST $ (`makeGradientBackground` small)

    -- Blend the background and the foreground.
    clouds <- makeClouds logger fp
    let image :: MIO.Image M.D RGB Float
        image = M.zipWith (\cloud bg -> interpolate cloud bg white) clouds gradientBackground
    MIO.writeImage "2d.png" $ M.map toWord8 image

main3d :: IO ()
main3d = do
    let numkeyframes =  32
        nonkeyframes =   2
        gifDelay     =   4

        height   = 200
        width    = 300
        small    = M.Ix3 numkeyframes height width
        world    = M.Sz (1 + 2 * small)
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
            , fpThreshold  = 0.6
            }

    -- Logging
    lock <- MVar.newMVar ()
    let logger = MVar.modifyMVar_ lock . const . IO.hPutStrLn IO.stderr

    cuby <- makeClouds logger fp

    logger "Generating keyframes..."
    let keyframes :: M.Array M.D M.Ix1 (M.Array M.M M.Ix2 Float)
        keyframes = M.map (cuby M.!>) (0 ..: numkeyframes)

        frames :: [JP.Image Word8]
        frames = M.toList $ M.map MIO.toJPImageY8 $
                 interpolateFrames nonkeyframes $ M.compute keyframes

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
                , JP.gfDelay       = gifDelay
                , JP.gfDisposal    = JP.DisposalDoNot
                , JP.gfPixels      = frame
                }
        }

main :: IO ()
main = main2d >> main3d
