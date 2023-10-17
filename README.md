# noise-ts

A simple library for 2D and 3D simplex noise and perlin noise in TypeScript.

This is a TypeScript implementation of [`josephg/noisejs`](https://github.com/josephg/noisejs). Slightly modified from [`alanko0511/noisejs-typescript`](https://github.com/alanko0511/noisejs-typescript).

## Installation
```python
npm i noise-ts
```

## How to use
```js
import Noise from 'noise-ts';

// Seed value is optional, default is 0.
const seed = Math.random();
const noise = new Noise(seed);

// Call the noise functions to get the noise value for that coordinates.
noise.simplex2(x, y);
noise.simplex3(x, y, z);

noise.perlin2(x, y);
noise.perlin3(x, y, z);
```
Please read the original repo's README for more information.
