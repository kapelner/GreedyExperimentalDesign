# Install Notes

## Configure Script Permissions

This package includes `configure` and `configure.win` scripts to remove stale
compiled artifacts before building. These scripts must be executable on Unix
systems. If you build from a zip or copy that drops file permissions, fix it with:

```sh
chmod +x configure configure.win
```

Then rerun:

```sh
R CMD INSTALL GreedyExperimentalDesign
```

## Windows Cleanup

On Windows, `configure.win` is a no-op to avoid install-time hangs. If you see
architecture mismatch errors (e.g., "not a valid Win32 application"), manually
remove stale artifacts in `src/` before reinstalling:

```sh
del /Q src\\*.o src\\*.so src\\*.dll src\\*.dylib
```
