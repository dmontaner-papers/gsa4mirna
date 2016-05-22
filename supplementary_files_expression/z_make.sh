#!/bin/bash

pandoc -s -S -f markdown -t html -o expressed_genes_report.html expressed_genes_report.md

pandoc -s -S -V papersize:a4paper -V geometry:"left=2cm, right=2cm, top=2cm, bottom=2cm" -f markdown -t latex -o expressed_genes_report.pdf expressed_genes_report.md

#pandoc -s -S -f markdown -t latex -o expressed_genes_report.tex expressed_genes_report.md
