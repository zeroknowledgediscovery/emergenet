#!/bin/bash

pdoc --html emergenet/ -o docs/ -c latex_math=True -f --template-dir docs/dark_templates

cp -r docs/emergenet/* docs
