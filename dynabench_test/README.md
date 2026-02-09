# SINDy on Dynabench Advection

## Usage

1. Download advection data:

    ```
    python -c 'from dynabench.dataset import download_equation; \
    download_equation("advection", structure="grid", resolution="medium")'
    ```

2. Run:

    ```uv run python3 validate_sindy.py```

## Notes
- Custom data loader written due to issues with `DynabenchIterator` and `DynabenchSimulationIterator`
- Only tested with low-resolution advection data on a grid. Unstructured (`cloud`) is not compatible with `PDELibrary`, and I could not get `WeakPDELibrary` to work (possible issues with source code).
- All animations/plots saved to `output/`
