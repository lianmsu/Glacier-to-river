#!/bin/bash
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------


minced -spacers ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/5455MAG_contig_withcas.fa ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/1minced_res/5455MAG_contig_withcas_minced.crisprs ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/1minced_res/5455MAG_contig_withcas_minced.gff

