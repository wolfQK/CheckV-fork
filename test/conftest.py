import shutil

import pytest


def pytest_addoption(parser):
    parser.addoption("--database", action="store", help="Path to CheckV's database")


@pytest.fixture
def database(request):
    return request.config.getoption("--database")


def pytest_sessionfinish(session, exitstatus):
    shutil.rmtree("test/output_files")
