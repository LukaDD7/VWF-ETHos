from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable


@dataclass(frozen=True)
class ContextualExcerpt:
    text: str
    before: str
    matched: str
    after: str
    start: int
    end: int
    nearest_variant_term: str | None
    nearest_variant_distance: int | None
    variant_linked: bool


def normalize_whitespace(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip()


def find_term_positions(text: str, term: str) -> list[tuple[int, int]]:
    if not term:
        return []
    return [
        (match.start(), match.end())
        for match in re.finditer(re.escape(term), text, flags=re.IGNORECASE)
    ]


def nearest_variant_position(
    text: str,
    start: int,
    end: int,
    variant_terms: Iterable[str],
) -> tuple[str | None, int | None, tuple[int, int] | None]:
    best_term: str | None = None
    best_distance: int | None = None
    best_position: tuple[int, int] | None = None
    for term in variant_terms:
        if not term:
            continue
        for variant_start, variant_end in find_term_positions(text, term):
            if variant_start <= end and start <= variant_end:
                distance = 0
            elif variant_start >= end:
                distance = variant_start - end
            else:
                distance = start - variant_end
            if best_distance is None or distance < best_distance:
                best_term = term
                best_distance = distance
                best_position = (variant_start, variant_end)
    return best_term, best_distance, best_position


def build_contextual_excerpt(
    *,
    text: str,
    start: int,
    end: int,
    term: str,
    variant_terms: Iterable[str],
    context_before_chars: int = 600,
    context_after_chars: int = 900,
    variant_link_radius: int = 1500,
) -> ContextualExcerpt:
    before_start = max(0, start - context_before_chars)
    after_end = min(len(text), end + context_after_chars)
    before = normalize_whitespace(text[before_start:start])
    matched = normalize_whitespace(text[start:end])
    after = normalize_whitespace(text[end:after_end])
    nearest_term, nearest_distance, _ = nearest_variant_position(text, start, end, variant_terms)
    variant_linked = nearest_distance is not None and nearest_distance <= variant_link_radius
    excerpt = " ".join(part for part in (before, matched, after) if part)
    if before_start > 0:
        excerpt = "… " + excerpt
    if after_end < len(text):
        excerpt += " …"
    return ContextualExcerpt(
        text=excerpt,
        before=before,
        matched=matched,
        after=after,
        start=start,
        end=end,
        nearest_variant_term=nearest_term,
        nearest_variant_distance=nearest_distance,
        variant_linked=variant_linked,
    )


def best_match_near_variant(
    *,
    text: str,
    term: str,
    variant_terms: Iterable[str],
) -> tuple[int, int] | None:
    positions = find_term_positions(text, term)
    if not positions:
        return None
    variant_terms = list(variant_terms)
    if not variant_terms:
        return positions[0]
    return min(
        positions,
        key=lambda position: nearest_variant_position(text, position[0], position[1], variant_terms)[1] or 10**9,
    )
