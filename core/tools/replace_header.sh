#!/bin/sh

infile=$1

perl ./replace_header.pl < $infile > tmp.file
mv tmp.file $infile
