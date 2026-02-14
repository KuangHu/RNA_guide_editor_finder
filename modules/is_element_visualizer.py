"""
IS Element Visualizer Module

Generates per-element PNG diagrams using dna_features_viewer showing IS element
structure: flanking regions, ORFs (directional arrows), noncoding regions, and
alignment hit positions.

Author: Kuang Hu
Date: 2026-02-10
"""

import logging
import os
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

logger = logging.getLogger(__name__)

# Color palette for IS families
IS_FAMILY_COLORS = {
    "IS110": "#e41a1c",
    "IS3": "#377eb8",
    "IS5": "#4daf4a",
    "IS21": "#984ea3",
    "IS630": "#ff7f00",
    "IS1182": "#a65628",
    "IS481": "#f781bf",
    "IS66": "#999999",
    "IS256": "#66c2a5",
    "ISL3": "#fc8d62",
    "IS30": "#8da0cb",
    "IS4": "#e78ac3",
    "IS1": "#a6d854",
    "IS1595": "#ffd92f",
    "IS91": "#e5c494",
    "IS701": "#b3b3b3",
    "IS1380": "#1b9e77",
    "IS982": "#d95f02",
    "IS6": "#7570b3",
    "IS1634": "#e7298a",
    "ISNCY": "#66a61e",
    "unclassified": "#cccccc",
}
DEFAULT_ORF_COLOR = "#bbbbbb"
NC_COLOR = "#e0e0e0"
FLANKING_UP_COLOR = "#aed6f1"     # light blue for upstream flanking
FLANKING_DOWN_COLOR = "#f9e79f"   # light yellow for downstream flanking
ALIGNMENT_UPSTREAM_COLOR = "#d62728"    # red for upstream flanking hits
ALIGNMENT_DOWNSTREAM_COLOR = "#1f77b4"  # blue for downstream flanking hits


def _deduplicate_alignments(alignments: List[Dict]) -> List[Dict]:
    """Remove duplicate alignments based on position, length, source, and orientation."""
    seen = set()
    deduped = []
    for aln in alignments:
        ungapped = aln.get("ungapped", {})
        key = (
            aln.get("flanking_source"),
            aln.get("noncoding_start"),
            ungapped.get("pos_in_noncoding"),
            ungapped.get("length"),
            ungapped.get("orientation"),
        )
        if key not in seen:
            seen.add(key)
            deduped.append(aln)
    return deduped


class ISElementVisualizer:
    """Generate diagrams of IS elements with flanking regions, ORFs, noncoding regions, and alignment hits."""

    def visualize_element(
        self,
        element: Dict,
        alignments: List[Dict],
        orf_annotations: Optional[List[Dict]] = None,
    ) -> plt.Figure:
        """Build a dna_features_viewer diagram for one IS element.

        The diagram spans from -upstream_length to is_length + downstream_length,
        showing flanking regions on either side of the IS element body.

        Args:
            element: dict from is_elements.json with is_element, orf_annotation, etc.
            alignments: list of passed alignment dicts for this element.
            orf_annotations: optional HMM annotation list to color ORFs by family.

        Returns:
            matplotlib Figure.
        """
        is_length = element.get("is_element", {}).get("length", 0)
        orf_ann = element.get("orf_annotation", {})
        orfs = orf_ann.get("orfs", [])
        nc_regions = orf_ann.get("noncoding_regions", [])

        # Flanking region lengths
        up_len = element.get("flanking_upstream", {}).get("length", 0)
        down_len = element.get("flanking_downstream", {}).get("length", 0)

        # Coordinate system: flanking_upstream occupies [-up_len, 0),
        # IS element occupies [0, is_length), downstream occupies [is_length, is_length+down_len)
        total_length = up_len + is_length + down_len
        offset = up_len  # shift all IS-internal coords by this amount

        # Build annotation lookup: orf_id -> is_family
        family_by_orf = {}
        if orf_annotations:
            for ann in orf_annotations:
                family_by_orf[ann["orf_id"]] = ann.get("is_family", "unclassified")

        features = []

        # 0. Flanking regions
        if up_len > 0:
            features.append(GraphicFeature(
                start=0,
                end=up_len,
                strand=0,
                color=FLANKING_UP_COLOR,
                label=f"upstream ({up_len}bp)",
                linewidth=0.5,
            ))
        if down_len > 0:
            features.append(GraphicFeature(
                start=offset + is_length,
                end=offset + is_length + down_len,
                strand=0,
                color=FLANKING_DOWN_COLOR,
                label=f"downstream ({down_len}bp)",
                linewidth=0.5,
            ))

        # 1. Noncoding regions — light gray rectangles
        for nc in nc_regions:
            nc_start = offset + nc["start"] - 1  # 1-based to 0-based + offset
            nc_end = offset + nc["end"]
            nc_type = nc.get("type", "")
            label = nc_type.replace("_", " ") if nc_type else "NC"
            nc_len = nc_end - nc_start
            features.append(GraphicFeature(
                start=nc_start,
                end=nc_end,
                strand=0,
                color=NC_COLOR,
                label=label if nc_len >= 30 else None,
                linewidth=0.5,
            ))

        # 2. ORFs — colored directional arrows
        for orf in orfs:
            orf_start = offset + orf["start"] - 1  # 1-based to 0-based + offset
            orf_end = offset + orf["end"]
            strand = +1 if orf.get("strand", "+") == "+" else -1
            orf_id = orf.get("orf_id", "")
            length_aa = orf.get("length_aa", 0)

            short_id = orf_id.rsplit("_", 1)[-1] if "_" in orf_id else orf_id
            label = f"orf{short_id} ({length_aa}aa)"

            family = family_by_orf.get(orf_id)
            color = IS_FAMILY_COLORS.get(family, DEFAULT_ORF_COLOR) if family else DEFAULT_ORF_COLOR

            features.append(GraphicFeature(
                start=orf_start,
                end=orf_end,
                strand=strand,
                color=color,
                label=label,
                linewidth=1,
            ))

        # 3. Alignment hits — deduplicated, colored markers on both flanking and noncoding sides
        deduped = _deduplicate_alignments(alignments)
        alignment_pairs = []  # collect (nc_midpoint, flank_midpoint, color) for connecting lines
        for aln in deduped:
            ungapped = aln.get("ungapped", {})
            nc_start_1based = aln.get("noncoding_start", 0)
            pos_in_nc = ungapped.get("pos_in_noncoding", 0)
            pos_in_flank = ungapped.get("pos_in_flanking", 0)
            aln_len = ungapped.get("length", 0)

            flanking = aln.get("flanking_source", "?")
            color = ALIGNMENT_UPSTREAM_COLOR if flanking == "upstream" else ALIGNMENT_DOWNSTREAM_COLOR
            orientation = ungapped.get("orientation", "forward")
            arrow = "\u2191" if flanking == "upstream" else "\u2193"
            ori_tag = "" if orientation == "forward" else " rc"
            label = f"{arrow}{aln_len}bp{ori_tag}"

            # Noncoding-side position (0-based diagram coords)
            nc_hit_start = offset + (nc_start_1based - 1) + pos_in_nc
            nc_hit_end = nc_hit_start + aln_len

            features.append(GraphicFeature(
                start=nc_hit_start,
                end=nc_hit_end,
                strand=0,
                color=color,
                label=label,
                linewidth=1.5,
                linecolor=color,
            ))

            # Flanking-side position (0-based diagram coords)
            if flanking == "upstream" and up_len > 0:
                flank_hit_start = pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            elif flanking == "downstream" and down_len > 0:
                flank_hit_start = offset + is_length + pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            else:
                continue  # no flanking region to draw in

            features.append(GraphicFeature(
                start=flank_hit_start,
                end=flank_hit_end,
                strand=0,
                color=color,
                label=None,
                linewidth=1.5,
                linecolor=color,
            ))

            # Save midpoints for connecting lines
            nc_mid = (nc_hit_start + nc_hit_end) / 2
            flank_mid = (flank_hit_start + flank_hit_end) / 2
            alignment_pairs.append((nc_mid, flank_mid, color))

        record = GraphicRecord(sequence_length=total_length, features=features)
        ax, _ = record.plot(figure_width=max(8, total_length / 150))

        # Draw connecting lines between flanking-side and noncoding-side hits
        for nc_mid, flank_mid, color in alignment_pairs:
            ax.annotate(
                "", xy=(nc_mid, -0.4), xytext=(flank_mid, -0.4),
                arrowprops=dict(
                    arrowstyle="-",
                    color=color,
                    alpha=0.35,
                    linewidth=1.5,
                    connectionstyle="arc3,rad=0.3",
                ),
            )

        # Draw IS element boundary lines
        ax.axvline(x=offset, color="black", linewidth=1.5, linestyle="--", alpha=0.6)
        ax.axvline(x=offset + is_length, color="black", linewidth=1.5, linestyle="--", alpha=0.6)

        fig = ax.figure
        return fig

    def save_element_png(
        self,
        element: Dict,
        alignments: List[Dict],
        output_path: str,
        dpi: int = 150,
        figure_width: int = 12,
        orf_annotations: Optional[List[Dict]] = None,
    ):
        """Generate and save a PNG diagram for one IS element.

        Args:
            element: dict from is_elements.json.
            alignments: passed alignment dicts for this element.
            output_path: path to save the PNG.
            dpi: resolution.
            figure_width: figure width in inches.
            orf_annotations: optional HMM annotations for coloring ORFs.
        """
        fig = self.visualize_element(element, alignments, orf_annotations)

        # Build title
        is_id = element.get("is_id", "unknown")
        is_length = element.get("is_element", {}).get("length", 0)

        family = "unknown"
        if orf_annotations:
            families = {a.get("is_family", "") for a in orf_annotations if a.get("is_family")}
            families.discard("")
            if families:
                family = ", ".join(sorted(families))

        deduped_count = len(_deduplicate_alignments(alignments))
        title = f"{is_id}  ({is_length} bp)"
        if family != "unknown":
            title += f"  [{family}]"
        title += f"  \u2014  {deduped_count} alignment(s)"

        fig.axes[0].set_title(title, fontsize=10, pad=10)
        fig.set_size_inches(figure_width, fig.get_size_inches()[1])

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def batch_visualize(
        self,
        elements_with_alignments: List[Dict],
        output_dir: str,
        dpi: int = 150,
        figure_width: int = 12,
    ):
        """Generate PNGs for multiple elements.

        Args:
            elements_with_alignments: list of dicts, each with keys:
                "element" (from is_elements.json),
                "alignments" (list of passed alignment dicts),
                "orf_annotations" (optional list of HMM annotations).
            output_dir: directory to write PNGs.
            dpi: resolution.
            figure_width: figure width in inches.
        """
        os.makedirs(output_dir, exist_ok=True)
        total = len(elements_with_alignments)

        for i, item in enumerate(elements_with_alignments, 1):
            element = item["element"]
            alignments = item["alignments"]
            orf_annotations = item.get("orf_annotations")
            is_id = element.get("is_id", f"unknown_{i}")

            output_path = os.path.join(output_dir, f"{is_id}.png")
            try:
                self.save_element_png(
                    element, alignments, output_path,
                    dpi=dpi, figure_width=figure_width,
                    orf_annotations=orf_annotations,
                )
                if i % 50 == 0 or i == total:
                    logger.info("Progress: %d/%d elements", i, total)
            except Exception:
                logger.exception("Failed to visualize %s", is_id)
