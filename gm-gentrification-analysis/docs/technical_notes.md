
# Technical Notes

- CRS: All spatial distance/area operations use EPSG:27700 (British National Grid).
- Validity: Use `st_make_valid` (R) / overlay with clean geometries (Python) before intersections.
- Performance: Enable spatial indexing; process large layers in chunks; cache intermediates.
- Reproducibility: See `requirements.txt`, `R-packages.txt`, `Stata-packages.txt`.
