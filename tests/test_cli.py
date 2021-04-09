#!/usr/bin/env python
import pytest


def test_help(script_runner):
    ret = script_runner.run("aMGSIM", "ancient-genomes", "-h")
    assert ret.success
