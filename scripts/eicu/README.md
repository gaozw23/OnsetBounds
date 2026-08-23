# eICU pipeline

This directory contains portable interfaces for eICU v2.0. It intentionally excludes raw hospital identifiers, patient records, fitted objects, preprocessing arrays, thresholds, and bootstrap draws. Authorized local data paths and TRAIN-derived constants must be supplied outside the repository.

The primary inference scheme resamples hospitals and then patients within sampled hospitals. The conditional-on-observed-hospitals analysis is documented as a sensitivity analysis and is not invoked by public CI.
