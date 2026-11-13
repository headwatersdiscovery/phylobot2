from django import template

register = template.Library()


@register.filter
def duration(total_seconds: float | None) -> str:
    """Seconds as '45 s', '3 min 5 s' or '2 h 10 min'; '–' if unknown."""
    if total_seconds is None:
        return "–"
    seconds = int(total_seconds)
    if seconds < 60:
        return f"{seconds} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes} min {seconds} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


@register.filter
def step_duration(total_seconds: float | None) -> str:
    """As `duration`, but '<1 s' for the quickest pipeline steps; '' if unknown."""
    if total_seconds is None:
        return ""
    return "<1 s" if total_seconds < 1 else duration(total_seconds)


@register.filter
def aligner_of(key: str) -> str:
    """The alignment a result file belongs to: its path's first part (mafft/...)."""
    return key.split("/")[0]
