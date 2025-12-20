# Makefile for reproducible grading tasks
.PHONY: benchmarks report verify

benchmarks:
	@echo "Running full benchmarks + plots..."
	python3 scripts/benchmark_and_plot.py

report:
	@echo "Generating concise report (markdown)..."
	python3 scripts/generate_report_md.py

verify:
	@echo "Verifying key outputs..."
	@test -f results/benchmark_report.md && echo "OK: results/benchmark_report.md exists" || echo "MISSING: results/benchmark_report.md"
	# `benchmarkanalysis.md` is now included in `benchmark_report.md`
	@test -f results/fst_comparison.png && echo "OK: results/fst_comparison.png exists" || echo "MISSING: results/fst_comparison.png"
