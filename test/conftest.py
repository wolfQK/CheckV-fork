import shutil

import pytest


def pytest_addoption(parser):
    parser.addoption("--database", action="store", help="Path to CheckV's database")
    parser.addoption("--threads", default=1, type=int, action="store", help="Threads to use")


@pytest.fixture
def database(request):
    return request.config.getoption("--database")

@pytest.fixture
def threads(request):
    return request.config.getoption("--threads")

def pytest_sessionfinish(session, exitstatus):
    shutil.rmtree("test/output_files")
