from __future__ import annotations


def parse_candidate_text(text: str) -> tuple[str, ...]:
    values: list[str] = []
    normalized = str(text or "").replace("\uFF0C", ",").replace("\uFF1B", ",").replace(";", ",")
    for raw in normalized.splitlines():
        for item in raw.split(","):
            cleaned = item.strip()
            if cleaned and cleaned not in values:
                values.append(cleaned)
    return tuple(values)


def merge_candidate_text(existing: str, additions: list[str]) -> str:
    values = list(parse_candidate_text(existing))
    for item in additions:
        if item and item not in values:
            values.append(item)
    return "\n".join(values)


__all__ = ["merge_candidate_text", "parse_candidate_text"]
