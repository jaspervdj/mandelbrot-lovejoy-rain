{-# LANGUAGE DeriveFunctor       #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import           Control.Monad                         (replicateM)
import qualified Data.Array.Accelerate                 as A
import qualified Data.Array.Accelerate.Data.Colour.RGB as A.RGB
import qualified Data.Array.Accelerate.IO.Codec.BMP    as A.BMP
import qualified Data.Array.Accelerate.LLVM.Native     as A.CPU
import qualified Test.QuickCheck                       as QC

data Pulse = Pulse
    { pCenter    :: (Float, Float)
    , pAmplitude :: Float
    , pRadius    :: Radius Float
    }

arbitraryPulse :: Float -> Float -> QC.Gen Pulse
arbitraryPulse alpha l = do
    center  <- (,) <$> QC.choose (0.0, l) <*> QC.choose (0.0, l)
    rhoInv <- QC.choose (0.0, 1.0)
    let rho   =  if rhoInv == 0.0 then 1.0 else 1.0 / (1.0 - rhoInv)
        radius = rho / 2.0
        amp    = rho ** (1.0 / alpha)
    ampSign <- QC.elements [-1.0, 1.0]
    return $ Pulse center (ampSign * amp) (Radius radius)

newtype Radius   f = Radius   f deriving (Functor)
newtype Distance f = Distance f deriving (Functor)

data PulseShape f
    = Rectangular
    | Smooth f    -- s
    | Annuli f f  -- lambda, s

pulseAt
    :: PulseShape (A.Exp Float) -> Distance (A.Exp Float) -> Radius (A.Exp Float)
    -> A.Exp Float
pulseAt Rectangular       (Distance u) (Radius r) =
    A.max (0.0 :: A.Exp Float) $
    A.fromIntegral (A.ceiling (1.0 - u / r) :: A.Exp Int)
pulseAt (Smooth s)        (Distance u) (Radius r) = exp (-(u / r) ** (2 * s))
pulseAt (Annuli lambda s) (Distance u) (Radius r) =
    let lambda' = sqrt (lambda * lambda - 1.0)
        delta   = (lambda + lambda') / 2.0
        sigma   = (lambda - lambda') / 2.0 in
    exp (-(((u*u) / (r*r) - (delta*delta)) / (sigma*sigma)) ** (2 * s))

generatePulse
    :: (Int, Int)
    -> PulseShape (A.Exp Float)
    -> Pulse
    -> A.Acc (A.Array A.DIM2 Float)
generatePulse (w, h) shape pulse = A.compute $ A.generate
    (A.constant (A.Z A.:. w A.:. h))
    (\idx ->
        let (A.Z A.:. ix A.:. iy) = A.unlift idx
            x                     = A.fromIntegral ix
            y                     = A.fromIntegral iy
            cx                    = A.constant (fst $ pCenter pulse)
            cy                    = A.constant (snd $ pCenter pulse)
            dx                    = x - cx
            dy                    = y - cy
            dist                  = Distance (sqrt $ dx * dx + dy * dy) in
        pulseAt shape dist (A.constant <$> pRadius pulse))

doubleToColour :: A.Exp Float -> A.Exp A.RGB.Colour
doubleToColour x = A.RGB.rgb x x x

normalize
    :: A.Acc (A.Array A.DIM2 Float) -> A.Acc (A.Array A.DIM2 Float)
normalize arr =
    let flat = A.flatten arr
        mini = A.replicate (A.shape arr) (A.minimum flat)
        maxi = A.replicate (A.shape arr) (A.maximum flat) in
    A.zipWith3 (\x mi ma -> (x - mi) / (ma - mi)) arr mini maxi

main :: IO ()
main = do
    pulses <- QC.generate $ replicateM 1600 $ arbitraryPulse (5 / 3) (400)
    let accs  = map (generatePulse (400, 400) Rectangular) pulses
        summy = normalize $ foldl1 (A.zipWith (+)) accs
    A.BMP.writeImageToBMP "test.bmp" $ A.CPU.run $
        A.map A.RGB.packRGB $
        A.map doubleToColour $
        summy
