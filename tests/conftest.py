"""Shared pytest configuration for application tests."""


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: exercises the installed or command-line application")
    config.addinivalue_line("markers", "slow: long-running application regression test")
