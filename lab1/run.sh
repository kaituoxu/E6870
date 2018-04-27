#!/bin/bash

./lab1 --audio_file p1in.dat --feat_file temp \
  --frontend.window true --window.hamming true \
  --frontend.fft false \
  --frontend.melbin false \
  --frontend.dct false
