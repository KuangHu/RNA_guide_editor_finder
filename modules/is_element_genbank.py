"""
IS Element GenBank Export Module

Generates annotated GenBank (.gbk) files for IS elements with flanking regions,
ORFs, noncoding regions, and alignment hit positions. Output files are loadable
in SnapGene, Benchling, and other sequence viewers.

Author: Kuang Hu
Date: 2026-02-13
"""

import logging
import os
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .is_element_visualizer import _deduplicate_alignments

logger = logging.getLogger(__name__)


class ISElementGenBank:
    """Generate annotated GenBank files for IS elements."""

    def build_record(
        self,
        element: Dict,
        alignments: List[Dict],
        orf_annotations: Optional[List[Dict]] = None,
    ) -> SeqRecord:
        """Build a BioPython SeqRecord for one IS element with annotated features.

        Coordinate system matches the visualizer: upstream flanking first,
        then IS body, then downstream flanking.

        Args:
            element: dict from is_elements.json with is_element, orf_annotation, etc.
            alignments: list of passed alignment dicts for this element.
            orf_annotations: optional HMM annotation list to annotate ORFs by family.

        Returns:
            Bio.SeqRecord.SeqRecord with features.
        """
        is_id = element.get("is_id", "unknown")
        is_data = element.get("is_element", {})
        is_seq = is_data.get("sequence", "")
        is_length = is_data.get("length", len(is_seq))

        up_seq = element.get("flanking_upstream", {}).get("sequence", "")
        down_seq = element.get("flanking_downstream", {}).get("sequence", "")
        up_len = len(up_seq)
        down_len = len(down_seq)

        full_seq = up_seq + is_seq + down_seq

        orf_ann = element.get("orf_annotation", {})
        orfs = orf_ann.get("orfs", [])
        nc_regions = orf_ann.get("noncoding_regions", [])

        # Build ORF annotation lookup: orf_id -> annotation dict
        ann_by_orf = {}
        if orf_annotations:
            for ann in orf_annotations:
                ann_by_orf[ann["orf_id"]] = ann

        # Determine IS family for description
        family = ""
        if orf_annotations:
            families = {a.get("is_family", "") for a in orf_annotations if a.get("is_family")}
            families.discard("")
            if families:
                family = ", ".join(sorted(families))

        features = []

        # 1. Flanking regions — misc_feature
        if up_len > 0:
            features.append(SeqFeature(
                FeatureLocation(0, up_len, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": ["upstream_flanking"],
                    "note": [f"{up_len}bp upstream flanking region"],
                },
            ))
        if down_len > 0:
            features.append(SeqFeature(
                FeatureLocation(up_len + is_length, up_len + is_length + down_len, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": ["downstream_flanking"],
                    "note": [f"{down_len}bp downstream flanking region"],
                },
            ))

        # 2. IS element boundary — misc_feature
        is_note = f"{is_length}bp IS element"
        if family:
            is_note += f" [{family}]"
        sample = element.get("sample", "")
        if sample:
            is_note += f" sample={sample}"
        features.append(SeqFeature(
            FeatureLocation(up_len, up_len + is_length, strand=0),
            type="misc_feature",
            qualifiers={
                "label": [is_id],
                "note": [is_note],
            },
        ))

        # 3. ORFs — CDS
        for orf in orfs:
            orf_start = up_len + orf["start"] - 1  # 1-based to 0-based + offset
            orf_end = up_len + orf["end"]
            strand = +1 if orf.get("strand", "+") == "+" else -1
            orf_id = orf.get("orf_id", "")
            length_aa = orf.get("length_aa", 0)
            protein_seq = orf.get("protein_sequence", "")

            qualifiers = {
                "label": [orf_id],
                "note": [f"{length_aa}aa"],
            }
            if protein_seq:
                qualifiers["translation"] = [protein_seq]

            # Add HMM annotation info if available
            ann = ann_by_orf.get(orf_id)
            if ann:
                ann_family = ann.get("is_family", "")
                e_value = ann.get("e_value")
                if ann_family:
                    note = qualifiers["note"][0]
                    note += f"; IS family: {ann_family}"
                    if e_value is not None:
                        note += f"; e-value: {e_value}"
                    qualifiers["note"] = [note]

            features.append(SeqFeature(
                FeatureLocation(orf_start, orf_end, strand=strand),
                type="CDS",
                qualifiers=qualifiers,
            ))

        # 4. Noncoding regions — misc_feature
        for nc in nc_regions:
            nc_start = up_len + nc["start"] - 1  # 1-based to 0-based + offset
            nc_end = up_len + nc["end"]
            nc_type = nc.get("type", "noncoding")
            features.append(SeqFeature(
                FeatureLocation(nc_start, nc_end, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": [nc_type],
                    "note": [f"{nc_end - nc_start}bp noncoding region"],
                },
            ))

        # 5. Alignment hits — misc_binding for both noncoding and flanking sides
        deduped = _deduplicate_alignments(alignments)
        for aln in deduped:
            ungapped = aln.get("ungapped", {})
            nc_start_1based = aln.get("noncoding_start", 0)
            pos_in_nc = ungapped.get("pos_in_noncoding", 0)
            pos_in_flank = ungapped.get("pos_in_flanking", 0)
            aln_len = ungapped.get("length", 0)

            flanking = aln.get("flanking_source", "?")
            orientation = ungapped.get("orientation", "forward")
            e_value = aln.get("e_value")
            confidence = aln.get("confidence")

            # Build note
            note_parts = [f"orientation={orientation}", f"flanking_source={flanking}"]
            if e_value is not None:
                note_parts.append(f"e_value={e_value}")
            if confidence:
                note_parts.append(f"confidence={confidence}")
            note = "; ".join(note_parts)

            label = f"aln_{flanking}_{aln_len}bp"

            # Noncoding-side hit
            nc_hit_start = up_len + (nc_start_1based - 1) + pos_in_nc
            nc_hit_end = nc_hit_start + aln_len
            features.append(SeqFeature(
                FeatureLocation(nc_hit_start, nc_hit_end, strand=0),
                type="misc_binding",
                qualifiers={
                    "label": [label + "_nc"],
                    "note": [note],
                },
            ))

            # Flanking-side hit
            if flanking == "upstream" and up_len > 0:
                flank_hit_start = pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            elif flanking == "downstream" and down_len > 0:
                flank_hit_start = up_len + is_length + pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            else:
                continue

            features.append(SeqFeature(
                FeatureLocation(flank_hit_start, flank_hit_end, strand=0),
                type="misc_binding",
                qualifiers={
                    "label": [label + "_flank"],
                    "note": [note],
                },
            ))

        # Build record
        description = f"{is_id} {is_length}bp IS element"
        if family:
            description += f" [{family}]"

        record = SeqRecord(
            Seq(full_seq),
            id=is_id,
            name=is_id[:16],
            description=description,
            features=features,
            annotations={
                "molecule_type": "DNA",
                "topology": "linear",
            },
        )

        return record

    def save_genbank(
        self,
        element: Dict,
        alignments: List[Dict],
        output_path: str,
        orf_annotations: Optional[List[Dict]] = None,
    ):
        """Generate and save a GenBank file for one IS element.

        Args:
            element: dict from is_elements.json.
            alignments: passed alignment dicts for this element.
            output_path: path to save the .gbk file.
            orf_annotations: optional HMM annotations for ORFs.
        """
        record = self.build_record(element, alignments, orf_annotations)
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        SeqIO.write(record, output_path, "genbank")

    def batch_export(
        self,
        elements_with_alignments: List[Dict],
        output_dir: str,
    ):
        """Generate GenBank files for multiple elements.

        Args:
            elements_with_alignments: list of dicts, each with keys:
                "element" (from is_elements.json),
                "alignments" (list of passed alignment dicts),
                "orf_annotations" (optional list of HMM annotations).
            output_dir: directory to write .gbk files.
        """
        os.makedirs(output_dir, exist_ok=True)
        total = len(elements_with_alignments)

        for i, item in enumerate(elements_with_alignments, 1):
            element = item["element"]
            alignments = item["alignments"]
            orf_annotations = item.get("orf_annotations")
            is_id = element.get("is_id", f"unknown_{i}")

            output_path = os.path.join(output_dir, f"{is_id}.gbk")
            try:
                self.save_genbank(
                    element, alignments, output_path,
                    orf_annotations=orf_annotations,
                )
                if i % 50 == 0 or i == total:
                    logger.info("GenBank progress: %d/%d elements", i, total)
            except Exception:
                logger.exception("Failed to export GenBank for %s", is_id)
