#! /bin/bash

emu abundance S1.fq --db /home/marshaag/emu_gitlab/emu_database --threads 16 --keep-counts --output-dir emu_results
emu abundance S2.fq --db /home/marshaag/emu_gitlab/emu_database --threads 16 --keep-counts --output-dir emu_results
emu abundance S3.fq --db /home/marshaag/emu_gitlab/emu_database --threads 16 --keep-counts --output-dir emu_results
emu abundance S4.fq --db /home/marshaag/emu_gitlab/emu_database --threads 16 --keep-counts --output-dir emu_results
emu abundance S5.fq --db /home/marshaag/emu_gitlab/emu_database --threads 16 --keep-counts --output-dir emu_results