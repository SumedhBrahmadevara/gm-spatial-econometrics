
def test_import():
    import importlib
    m = importlib.import_module("code.gm_tools_improved")
    assert hasattr(m, "compute_spatial_controls")
