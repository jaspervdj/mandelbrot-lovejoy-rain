{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE ViewPatterns  #-}
module Main where

import qualified Data.Array.Accelerate                 as A
import qualified Data.Array.Accelerate.Data.Colour.RGB as A.RGB
import qualified Data.Array.Accelerate.Data.Complex    as A
import qualified Data.Array.Accelerate.IO.Codec.BMP    as A.BMP
import qualified Data.Array.Accelerate.LLVM.Native     as A.CPU

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
    A.min (1.0 :: A.Exp Float) $ A.fromIntegral (A.floor (u / r) :: A.Exp Int)
pulseAt (Smooth s)        (Distance u) (Radius r) = exp (-(u / r) ** (2 * s))
pulseAt (Annuli lambda s) (Distance u) (Radius r) =
    let lambda' = sqrt (lambda * lambda - 1.0)
        delta   = (lambda + lambda') / 2.0
        sigma   = (lambda - lambda') / 2.0 in
    exp (-(((u*u) / (r*r) - (delta*delta)) / (sigma*sigma)) ** (2 * s))

data Pulse = Pulse
    { pCenter    :: (Float, Float)
    , pAmplitude :: Float
    , pRadius    :: Radius Float
    }

generatePulse
    :: (Int, Int)
    -> PulseShape (A.Exp Float)
    -> Pulse
    -> A.Acc (A.Array A.DIM2 Float)
generatePulse (w, h) shape pulse = A.generate
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
  where
    indexToTuple :: A.Exp A.DIM2 -> A.Exp (Float, Float)
    indexToTuple (A.unlift -> A.Z A.:. y A.:. x) =
        A.lift (A.fromIntegral y, A.fromIntegral x)

doubleToColour :: A.Exp Float -> A.Exp A.RGB.Colour
doubleToColour x = A.RGB.rgb x x x

main :: IO ()
main =
    A.BMP.writeImageToBMP "test.bmp" $ A.CPU.run $
    A.map A.RGB.packRGB $
    A.map doubleToColour $
    generatePulse (200, 200) (Annuli 1.01 1) (Pulse (100, 100) 1.0 (Radius 50))
