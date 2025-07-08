import pytest


def pytest_configure(config):
    configure_is_imported = True
    # filter expected warnings
    config.addinivalue_line(
        "filterwarnings",
        "ignore:No geometry type found, no coordinates will be added:UserWarning",
    )
    config.addinivalue_line(
        "filterwarnings",
        "ignore:deallocating CachingFileManager.*, but file is not already closed. "
        "This may indicate a bug.:RuntimeWarning",
    )
    config.addinivalue_line(
        "filterwarnings",
        "ignore:While building x, y, z coordinate arrays, an exception occured"
        ":boututils.boutwarnings.AlwaysWarning",
    )
    config.addinivalue_line(
        "filterwarnings",
        "ignore:Haven't decided how to write options file back out yet:UserWarning",
    )

    # register additional markers
    config.addinivalue_line(
        "markers",
        "long: long test, or one of many permutations (disabled by default)",
    )
    config.addinivalue_line(
        "markers",
        "flaky: xarray uses this mark, adding to avoid warnings about it not being "
        "defined",
    )
    config.addinivalue_line(
        "markers",
        "network: xarray uses this mark, adding to avoid warnings about it not "
        "being defined",
    )
    config.addinivalue_line(
        "markers",
        "slow: xarray uses this mark, adding to avoid warnings about it not being "
        "defined",
    )


def pytest_addoption(parser):
    # Add command line option '--long' for pytest, to be used to enable long tests
    addoption_is_imported = True
    parser.addoption(
        "--long",
        action="store_true",
        default=False,
        help="enable tests marked as 'long'",
    )


def pytest_collection_modifyitems(config, items):
    collection_is_imported = True
    if not config.getoption("--long"):
        # --long not given in cli: skip long tests
        print("\n    skipping long tests, pass '--long' to enable")
        skip_long = pytest.mark.skip(reason="need --long option to run")
        for item in items:
            if "long" in item.keywords:
                item.add_marker(skip_long)
