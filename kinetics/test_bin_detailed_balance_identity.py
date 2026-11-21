"""Synthetic algebraic acceptance control, NOT a physical 5D scattering calculation."""
import math
import random

def test_pair_identity():
    rng = random.Random(47)
    for _ in range(10000):
        ei, ej = sorted((rng.uniform(.1, 4), rng.uniform(.1, 4)))
        t = rng.uniform(.2, 3)
        fi, fj = (rng.uniform(1e-5, 1-1e-5) for _ in range(2))
        a = math.exp((ei-ej)/(2*t))
        x = a*fi*(1-fj)
        y = fj*(1-fi)/a
        sigma = (x-y)*math.log(x/y)
        assert sigma >= -1e-14
        mu = rng.uniform(-2, 2)
        fi_eq = 1/(1+math.exp((ei-mu)/t))
        fj_eq = 1/(1+math.exp((ej-mu)/t))
        forward = a*fi_eq*(1-fj_eq)
        reverse = fj_eq*(1-fi_eq)/a
        assert abs(forward-reverse) <= 2e-14
    print("10000 seeded synthetic detailed-balance/entropy pairs passed")

if __name__ == "__main__":
    test_pair_identity()
