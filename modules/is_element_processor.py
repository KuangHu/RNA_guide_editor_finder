"""
IS Element Processor Module

Adapter that reads IS element JSON data (flanking sequences + noncoding regions),
runs ShortAlignmentFinder between every flanking x noncoding combination, tracks
provenance, and maps output to FilterEngine format for downstream filtering.

Author: Kuang Hu
Date: 2026-02-09
"""

import json
import logging
from typing import Dict, List, Optional

from .short_alignment_finder import ShortAlignmentFinder
from .alignment_filter import FilterEngine

logger = logging.getLogger(__name__)


class ISElementProcessor:
    """Process IS element JSON data to find flanking-noncoding alignments.

    Reads IS element records, runs find_alignments_between for every
    flanking (upstream/downstream) x noncoding region combination, annotates
    results with provenance, and optionally filters through a FilterEngine.
    """

    def __init__(self,
                 min_length: int = 8,
                 max_mismatches: int = 0,
                 check_forward: bool = True,
                 check_revcomp: bool = True,
                 extend_with_gaps: bool = False,
                 filter_pipeline=None):
        """
        Parameters:
            min_length: Minimum alignment length for ShortAlignmentFinder.
            max_mismatches: Maximum mismatches allowed.
            check_forward: Check forward (direct) matches.
            check_revcomp: Check reverse complement matches.
            extend_with_gaps: Whether to extend ungapped hits with gapped alignment.
            filter_pipeline: Optional FilterPipeline for downstream filtering.
        """
        self.finder = ShortAlignmentFinder(
            min_length=min_length,
            max_mismatches=max_mismatches,
            check_forward=check_forward,
            check_revcomp=check_revcomp,
        )
        self.extend_with_gaps = extend_with_gaps
        self.filter_engine = FilterEngine(filter_pipeline) if filter_pipeline else None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def process_element(self, element: Dict) -> Dict:
        """Process one IS element: run all flanking x noncoding alignments.

        Returns a per-element result dict with alignments and optional filter results.
        """
        is_id = element.get("is_id", "unknown")
        sample = element.get("sample", "unknown")
        is_length = element.get("is_element", {}).get("length", 0)

        flanking_seqs, nc_regions = self._extract_sequences(element)

        all_alignments: List[Dict] = []

        for flanking_source, flanking_seq in flanking_seqs.items():
            if not flanking_seq:
                logger.warning("%s: empty %s flanking sequence, skipping", is_id, flanking_source)
                continue
            for nc_index, nc_region in enumerate(nc_regions):
                nc_seq = nc_region.get("sequence", "")
                if not nc_seq:
                    continue
                hits = self._run_alignments(
                    flanking_seq, nc_seq,
                    flanking_source, nc_region, nc_index, element,
                )
                all_alignments.extend(hits)

        result: Dict = {
            "is_id": is_id,
            "sample": sample,
            "is_element_length": is_length,
            "num_noncoding_regions": len(nc_regions),
            "total_hits": len(all_alignments),
            "alignments": all_alignments,
        }

        if self.filter_engine is not None and all_alignments:
            result["filtered"] = self._run_filter(all_alignments, element)

        return result

    def process_file(self, json_path: str) -> List[Dict]:
        """Load an is_elements.json file and process every element."""
        elements = load_is_elements(json_path)
        return [self.process_element(el) for el in elements]

    def process_batch(self, json_paths: List[str]) -> List[Dict]:
        """Process multiple is_elements.json files."""
        results: List[Dict] = []
        for path in json_paths:
            results.extend(self.process_file(path))
        return results

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _extract_sequences(self, element: Dict):
        """Pull flanking + noncoding sequences from an element dict.

        Returns:
            (flanking_seqs, nc_regions) where flanking_seqs is
            {"upstream": seq, "downstream": seq} and nc_regions is the
            list from orf_annotation.noncoding_regions.
        """
        flanking_seqs = {}
        for direction in ("upstream", "downstream"):
            key = f"flanking_{direction}"
            section = element.get(key, {})
            flanking_seqs[direction] = section.get("sequence", "")

        nc_regions = (
            element
            .get("orf_annotation", {})
            .get("noncoding_regions", [])
        )
        return flanking_seqs, nc_regions

    def _run_alignments(self, flanking_seq: str, noncoding_seq: str,
                        flanking_source: str, nc_region: Dict,
                        nc_index: int, element: Dict) -> List[Dict]:
        """Run find_alignments_between + optional gapped extension for one pair.

        Returns list of annotated alignment dicts.
        """
        hits = self.finder.find_alignments_between(flanking_seq, noncoding_seq)

        annotated: List[Dict] = []
        for hit in hits:
            entry = self._annotate_hit(
                hit, flanking_seq, noncoding_seq,
                flanking_source, nc_region, nc_index, element,
            )

            if self.extend_with_gaps:
                gapped = self.finder.extend_alignment_with_gaps(
                    flanking_seq, noncoding_seq, hit,
                )
                if gapped is not None:
                    entry["gapped"] = {
                        "pos_in_flanking": gapped["pos1"],
                        "pos_in_noncoding": gapped["pos2"],
                        "end_in_flanking": gapped["end1"],
                        "end_in_noncoding": gapped["end2"],
                        "length_in_flanking": gapped["length1"],
                        "length_in_noncoding": gapped["length2"],
                        "seq1_aligned": gapped["seq1_aligned"],
                        "seq2_aligned": gapped["seq2_aligned"],
                        "alignment_length": gapped["alignment_length"],
                        "matches": gapped["matches"],
                        "mismatches": gapped["mismatches"],
                        "gaps": gapped["gaps"],
                        "identity": gapped["identity"],
                        "alignment_string": gapped["alignment_string"],
                        "orientation": gapped["orientation"],
                    }
                else:
                    entry["gapped"] = None
            else:
                entry["gapped"] = None

            annotated.append(entry)

        return annotated

    def _annotate_hit(self, hit: Dict, flanking_seq: str, noncoding_seq: str,
                      flanking_source: str, nc_region: Dict,
                      nc_index: int, element: Dict) -> Dict:
        """Wrap a raw alignment hit with provenance metadata."""
        return {
            # Provenance
            "is_id": element.get("is_id", "unknown"),
            "sample": element.get("sample", "unknown"),

            # Flanking source
            "flanking_source": flanking_source,
            "flanking_length": len(flanking_seq),

            # Noncoding region tracking
            "noncoding_region_index": nc_index,
            "noncoding_region_type": nc_region.get("type", ""),
            "noncoding_start": nc_region.get("start", 0),
            "noncoding_end": nc_region.get("end", 0),
            "noncoding_length": nc_region.get("length", len(noncoding_seq)),

            # Ungapped alignment (always present)
            "ungapped": {
                "pos_in_flanking": hit["pos1"],
                "pos_in_noncoding": hit["pos2"],
                "length": hit["length"],
                "seq_flanking": hit["seq1"],
                "seq_noncoding": hit["seq2"],
                "mismatches": hit["mismatches"],
                "mismatch_positions": hit["mismatch_positions"],
                "orientation": hit["orientation"],
            },

            # Gapped alignment placeholder (filled later by caller)
            "gapped": None,
        }

    def _map_to_filter_format(self, alignment_result: Dict) -> Dict:
        """Convert an annotated alignment result to FilterEngine input format.

        The FilterEngine.filter_alignment expects a dict with at least:
          - aligned_sequence
          - non_coding_start / non_coding_end
          - mismatches / gaps / alignment_string
          - query_length / target_length
        """
        ungapped = alignment_result["ungapped"]
        gapped = alignment_result.get("gapped")

        if gapped is not None:
            aligned_seq = gapped["seq1_aligned"].replace("-", "")
            mismatches = gapped["mismatches"]
            gaps = gapped["gaps"]
            alignment_string = gapped["alignment_string"]
        else:
            aligned_seq = ungapped["seq_flanking"]
            mismatches = ungapped["mismatches"]
            gaps = 0
            alignment_string = ""

        return {
            "aligned_sequence": aligned_seq,
            "non_coding_start": alignment_result["noncoding_start"],
            "non_coding_end": alignment_result["noncoding_end"],
            "mismatches": mismatches,
            "gaps": gaps,
            "alignment_string": alignment_string,
            "query_length": alignment_result["flanking_length"],
            "target_length": alignment_result["noncoding_length"],
        }

    def _build_transposon_data(self, element: Dict) -> Dict:
        """Build transposon context dict for boundary artifact filter."""
        is_elem = element.get("is_element", {})
        return {
            "start": 1,
            "end": is_elem.get("length", 0),
        }

    def _run_filter(self, alignments: List[Dict], element: Dict) -> Dict:
        """Run FilterEngine on all alignments for one element."""
        transposon_data = self._build_transposon_data(element)

        total = len(alignments)
        passed_count = 0
        failed_count = 0
        by_confidence: Dict[str, int] = {"high": 0, "low": 0, "rejected": 0}
        results_list: List[Dict] = []

        for aln in alignments:
            filter_input = self._map_to_filter_format(aln)
            filter_result = self.filter_engine.filter_alignment(
                filter_input, transposon_data,
            )

            if filter_result["pass"]:
                passed_count += 1
                conf = filter_result.get("confidence", "high")
                by_confidence[conf] = by_confidence.get(conf, 0) + 1
            else:
                failed_count += 1
                by_confidence["rejected"] = by_confidence.get("rejected", 0) + 1

            results_list.append({
                "alignment": aln,
                "filter_result": filter_result,
            })

        return {
            "total": total,
            "passed": passed_count,
            "failed": failed_count,
            "by_confidence": by_confidence,
            "results": results_list,
        }


# ------------------------------------------------------------------
# Convenience functions
# ------------------------------------------------------------------

def load_is_elements(json_path: str) -> List[Dict]:
    """Load IS elements from a JSON file.

    The file should contain a JSON array of IS element records.
    """
    with open(json_path, "r") as fh:
        data = json.load(fh)
    if isinstance(data, list):
        return data
    raise ValueError(f"Expected JSON array in {json_path}, got {type(data).__name__}")


def process_is_element(element: Dict,
                       min_length: int = 8,
                       max_mismatches: int = 0,
                       check_forward: bool = True,
                       check_revcomp: bool = True,
                       extend_with_gaps: bool = False,
                       filter_pipeline=None) -> Dict:
    """Convenience: process a single IS element dict."""
    processor = ISElementProcessor(
        min_length=min_length,
        max_mismatches=max_mismatches,
        check_forward=check_forward,
        check_revcomp=check_revcomp,
        extend_with_gaps=extend_with_gaps,
        filter_pipeline=filter_pipeline,
    )
    return processor.process_element(element)


def process_is_elements_file(json_path: str,
                             min_length: int = 8,
                             max_mismatches: int = 0,
                             check_forward: bool = True,
                             check_revcomp: bool = True,
                             extend_with_gaps: bool = False,
                             filter_pipeline=None) -> List[Dict]:
    """Convenience: load and process an entire is_elements.json file."""
    processor = ISElementProcessor(
        min_length=min_length,
        max_mismatches=max_mismatches,
        check_forward=check_forward,
        check_revcomp=check_revcomp,
        extend_with_gaps=extend_with_gaps,
        filter_pipeline=filter_pipeline,
    )
    return processor.process_file(json_path)
