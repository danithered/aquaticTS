#!/bin/bash

cp src/abra.qmd $1
cd $1

quarto render abra.qmd -P dir:$(pwd)
rm abra.qmd
open ./abra.html

