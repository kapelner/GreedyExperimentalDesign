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
