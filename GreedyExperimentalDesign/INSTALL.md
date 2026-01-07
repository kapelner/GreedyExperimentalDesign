# Install Notes

## Manual Cleanup

If you see architecture mismatch errors (e.g., "not a valid Win32 application"),
remove stale compiled artifacts in `src/` before reinstalling.

macOS/Linux (from `GreedyExperimentalDesign/`):

```sh
rm -f src/*.o src/*.so src/*.dll src/*.dylib
```

Windows (from `GreedyExperimentalDesign\\`):

```sh
del /Q src\\*.o src\\*.so src\\*.dll src\\*.dylib
```
