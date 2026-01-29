import os
from pathlib import Path
import pytest

@pytest.fixture
def test_data_path():
    env_value = os.environ.get("ANGCAL_TEST_DATA")
    if not env_value:
        raise RuntimeError("Environment variable ANGCAL_TEST_DATA is not set or is empty")

    return Path(env_value)
