
pandoc -s -S -f markdown -t html -o godard_comparison_report.html godard_comparison_report.md

pandoc -s -S -f markdown -t latex -o godard_comparison_report.pdf godard_comparison_report.md
