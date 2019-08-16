#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""__init__ file for slap package, allowing import of
    LAEModel and Catalogue.

"""

from .model import LAEModel
from .catalogue import Catalogue

__all__ = ['LAEModel', 'Catalogue']
