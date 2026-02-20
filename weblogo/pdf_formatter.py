"""Native PDF renderer for WebLogo.

Generates valid PDF 1.4 directly from LogoData/LogoFormat,
with no external dependencies (no Ghostscript).

Written by Claude / Opus 4.6 (2026-02-16)
"""

from .color import Color
from .logo import LogoData, LogoFormat, std_units

# Exact character bounding boxes for Arial-BoldMT, extracted from PostScript
# via Ghostscript `charpath pathbbox` at 1000pt, normalized to fractions of em.
# Format: (lx, ly, ux, uy) — lower-left and upper-right of visual bbox.
_CHAR_BBOX = {
    "A": (0.000, 0.000, 0.712, 0.710),
    "B": (0.072, 0.000, 0.668, 0.710),
    "C": (0.047, -0.012, 0.666, 0.723),
    "D": (0.072, 0.000, 0.667, 0.710),
    "E": (0.072, 0.000, 0.612, 0.710),
    "F": (0.072, 0.000, 0.560, 0.710),
    "G": (0.048, -0.013, 0.711, 0.723),
    "H": (0.072, 0.000, 0.640, 0.710),
    "I": (0.067, 0.000, 0.211, 0.710),
    "J": (0.017, -0.013, 0.471, 0.710),
    "K": (0.074, 0.000, 0.715, 0.710),
    "L": (0.076, 0.000, 0.576, 0.704),
    "M": (0.070, 0.000, 0.755, 0.710),
    "N": (0.072, 0.000, 0.637, 0.710),
    "O": (0.043, -0.013, 0.732, 0.723),
    "P": (0.072, 0.000, 0.616, 0.710),
    "Q": (0.043, -0.072, 0.758, 0.723),
    "R": (0.072, 0.000, 0.711, 0.710),
    "S": (0.035, -0.013, 0.612, 0.723),
    "T": (0.022, 0.000, 0.585, 0.710),
    "U": (0.072, -0.013, 0.637, 0.710),
    "V": (0.000, 0.000, 0.662, 0.710),
    "W": (0.003, 0.000, 0.933, 0.710),
    "X": (0.001, 0.000, 0.660, 0.710),
    "Y": (-0.001, 0.000, 0.663, 0.710),
    "Z": (0.010, 0.000, 0.587, 0.710),
}
_DEFAULT_BBOX = (0.0, 0.0, 0.60, 0.710)

# Times-Roman I bounding box (for the I replacement hack — serifed I)
_TIMES_I_BBOX = (0.018, 0.000, 0.315, 0.662)

# Helvetica character widths (fractions of em), from standard AFM data.
# Used for computing text widths for centering/right-aligning labels.
_HELVETICA_WIDTHS = {
    " ": 0.278, "!": 0.278, '"': 0.355, "#": 0.556, "$": 0.556,
    "%": 0.889, "&": 0.667, "'": 0.191, "(": 0.333, ")": 0.333,
    "*": 0.389, "+": 0.584, ",": 0.278, "-": 0.333, ".": 0.278,
    "/": 0.278, "0": 0.556, "1": 0.556, "2": 0.556, "3": 0.556,
    "4": 0.556, "5": 0.556, "6": 0.556, "7": 0.556, "8": 0.556,
    "9": 0.556, ":": 0.278, ";": 0.278, "<": 0.584, "=": 0.584,
    ">": 0.584, "?": 0.556, "@": 1.015, "A": 0.667, "B": 0.667,
    "C": 0.722, "D": 0.722, "E": 0.667, "F": 0.611, "G": 0.778,
    "H": 0.722, "I": 0.278, "J": 0.500, "K": 0.667, "L": 0.556,
    "M": 0.833, "N": 0.722, "O": 0.778, "P": 0.667, "Q": 0.778,
    "R": 0.722, "S": 0.667, "T": 0.611, "U": 0.722, "V": 0.667,
    "W": 0.944, "X": 0.667, "Y": 0.667, "Z": 0.611, "[": 0.278,
    "\\": 0.278, "]": 0.278, "^": 0.469, "_": 0.556, "`": 0.333,
    "a": 0.556, "b": 0.556, "c": 0.500, "d": 0.556, "e": 0.556,
    "f": 0.278, "g": 0.556, "h": 0.556, "i": 0.222, "j": 0.222,
    "k": 0.500, "l": 0.222, "m": 0.833, "n": 0.556, "o": 0.556,
    "p": 0.556, "q": 0.556, "r": 0.333, "s": 0.500, "t": 0.278,
    "u": 0.556, "v": 0.500, "w": 0.722, "x": 0.500, "y": 0.500,
    "z": 0.500, "{": 0.334, "|": 0.260, "}": 0.334, "~": 0.584,
    "\u2032": 0.333,  # prime symbol width (approximate)
}


def _pdf_escape(s: str) -> str:
    """Escape special characters for PDF text strings."""
    return s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _string_width(s: str, font_size: float) -> float:
    """Compute text width in points using Helvetica metrics."""
    w = sum(_HELVETICA_WIDTHS.get(c, 0.556) for c in s)
    return w * font_size


def native_pdf_formatter(logodata: LogoData, logoformat: LogoFormat) -> bytes:
    """Generate a logo in PDF format.

    Produces a valid PDF 1.4 with Standard 14 fonts — no external tools required.
    """
    assert logoformat.logo_width is not None
    assert logoformat.logo_height is not None

    logo_width = logoformat.logo_width
    logo_height = logoformat.logo_height
    margin = logoformat.logo_margin

    # Build the content stream
    stream_parts: list = []

    # --- Title ---
    if logoformat.show_title and logoformat.logo_title:
        _draw_title(stream_parts, logoformat, logo_width, logo_height)

    # --- Logo label (top left) ---
    if logoformat.logo_label:
        _draw_logo_label(stream_parts, logoformat, logo_height)

    # --- X-axis label (bottom center) ---
    if logoformat.show_xaxis_label and logoformat.xaxis_label:
        _draw_xaxis_label(stream_parts, logoformat, logo_width)

    # --- Fineprint (bottom right) ---
    if logoformat.show_fineprint and logoformat.fineprint:
        _draw_fineprint(stream_parts, logoformat, logo_width)

    # --- Lines ---
    assert logoformat.lines_per_logo is not None
    assert logoformat.line_height is not None
    assert logoformat.line_margin_left is not None
    assert logoformat.line_margin_right is not None
    assert logoformat.line_margin_bottom is not None
    assert logoformat.line_margin_top is not None
    assert logoformat.title_height is not None
    assert logoformat.xaxis_label_height is not None
    assert logoformat.logo_start is not None
    assert logoformat.first_index is not None
    assert logoformat.logo_end is not None
    assert logoformat.char_width is not None

    conv_factor = std_units[logoformat.unit_name]

    seq_from = logoformat.logo_start - logoformat.first_index
    seq_to = logoformat.logo_end - logoformat.first_index + 1

    line_index = 0
    stack_in_line = 0

    for seq_index in range(seq_from, seq_to):
        stack_index = seq_index - seq_from

        if stack_index != 0 and (stack_index % logoformat.stacks_per_line) == 0:
            line_index += 1
            stack_in_line = 0

        # PDF coordinate system: origin bottom-left, y-up.
        # Line bottom-left in logo coordinates:
        line_x = margin
        line_bottom = (
            margin
            + logoformat.xaxis_label_height
            + (logoformat.lines_per_logo - 1 - line_index) * logoformat.line_height
        )

        # Content area within the line
        content_x = line_x + logoformat.line_margin_left
        content_y = line_bottom + logoformat.line_margin_bottom

        # Draw y-axis (once per line)
        if stack_in_line == 0:
            if logoformat.show_yaxis:
                _draw_yaxis(stream_parts, logoformat, content_x, content_y)

            if logoformat.show_xaxis and logoformat.show_ends:
                _draw_left_end(stream_parts, logoformat, content_x, content_y)

        # Stack position
        stack_x = content_x + stack_in_line * logoformat.stack_width
        stack_y = content_y  # bottom of the stack

        # Draw x-axis tic and number
        if logoformat.show_xaxis:
            _draw_xaxis_tic(
                stream_parts, logoformat, stack_x, stack_y,
                logoformat.annotate[seq_index]
            )

        # Calculate stack height in data units
        if conv_factor:
            assert logodata.entropy is not None
            stack_height_units = (
                logodata.entropy[seq_index] * std_units[logoformat.unit_name]
            )
        else:
            stack_height_units = 1.0

        # Sort symbols by frequency
        assert logodata.alphabet is not None
        assert logodata.counts is not None
        s = list(zip(logodata.counts[seq_index], logodata.alphabet))
        s.sort(key=lambda x: x[1])
        s.reverse()
        s.sort(key=lambda x: x[0])

        if not logoformat.reverse_stacks:
            s.reverse()

        C = float(sum(logodata.counts[seq_index]))

        if C > 0.0:
            fraction_width = 1.0
            assert logoformat.scale_width is not None
            if logoformat.scale_width:
                assert logodata.weight is not None
                fraction_width = float(logodata.weight[seq_index])

            # Current y position (bottom of next symbol, building upward)
            y_cursor = stack_y

            for rank, (count, symbol) in enumerate(s):
                interval = count * stack_height_units / C
                slot_height = interval / logoformat.yaxis_scale * logoformat.stack_height
                char_height_pts = slot_height - 2 * logoformat.stack_margin

                if char_height_pts <= 0.01:
                    y_cursor += slot_height
                    continue

                target_width = fraction_width * logoformat.char_width
                if logoformat.show_boxes:
                    target_width *= logoformat.shrink_fraction
                    char_height_pts_draw = char_height_pts * logoformat.shrink_fraction
                else:
                    char_height_pts_draw = char_height_pts

                assert logoformat.color_scheme is not None
                color = logoformat.color_scheme.symbol_color(seq_index, symbol, rank)

                char_x = stack_x + logoformat.stack_margin
                char_x += (1 - fraction_width) * logoformat.char_width / 2

                char_y_bottom = y_cursor + logoformat.stack_margin

                if logoformat.show_boxes:
                    char_x += target_width * (1 - logoformat.shrink_fraction) / 2
                    char_y_bottom += (char_height_pts - char_height_pts_draw) / 2
                    target_width = target_width * logoformat.shrink_fraction

                _draw_symbol(
                    stream_parts,
                    symbol,
                    char_x,
                    char_y_bottom,
                    target_width,
                    char_height_pts_draw if logoformat.show_boxes else char_height_pts,
                    color,
                )

                y_cursor += slot_height

        # Draw error bar
        if logodata.entropy_interval is not None and conv_factor and C > 0.0:
            low, high = logodata.entropy_interval[seq_index]
            assert logodata.entropy is not None
            center = logodata.entropy[seq_index]
            low *= conv_factor
            high *= conv_factor
            center *= conv_factor

            if high > logoformat.yaxis_scale:
                high = logoformat.yaxis_scale

            down = center - low
            up = high - center
            _draw_errorbar(
                stream_parts, logoformat, stack_x, stack_y, center, down, up
            )

        # Draw right end after last stack in a line
        stack_in_line += 1
        is_last_stack = stack_index == seq_to - seq_from - 1
        is_line_end = stack_in_line == logoformat.stacks_per_line

        if (is_last_stack or is_line_end) and logoformat.show_xaxis and logoformat.show_ends:
            right_end_x = content_x + stack_in_line * logoformat.stack_width
            _draw_right_end(stream_parts, logoformat, right_end_x, content_y)

    content_stream = "\n".join(stream_parts)
    return _assemble_pdf(content_stream, logo_width, logo_height)


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _draw_title(parts: list, fmt: LogoFormat, logo_width: int,
                logo_height: int) -> None:
    """Draw the logo title, centered at top."""
    title_width = _string_width(fmt.logo_title, fmt.title_fontsize)
    x = (logo_width - title_width) / 2.0
    # In PDF coords (bottom-left origin), top of logo is at logo_height
    y = logo_height - fmt.logo_margin - fmt.title_fontsize * 0.75
    parts.append("BT")
    parts.append(f"/F2 {fmt.title_fontsize} Tf")
    parts.append(f"{x:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(fmt.logo_title)}) Tj")
    parts.append("ET")


def _draw_logo_label(parts: list, fmt: LogoFormat, logo_height: int) -> None:
    """Draw logo label at top-left."""
    x = fmt.logo_margin
    y = logo_height - fmt.logo_margin - fmt.title_fontsize * 0.75
    parts.append("BT")
    parts.append(f"/F2 {fmt.title_fontsize} Tf")
    parts.append(f"{x:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(fmt.logo_label)}) Tj")
    parts.append("ET")


def _draw_xaxis_label(parts: list, fmt: LogoFormat, logo_width: int) -> None:
    """Draw x-axis label, bottom center."""
    assert fmt.xaxis_label_height is not None
    label_width = _string_width(fmt.xaxis_label, fmt.fontsize)
    x = (logo_width - label_width) / 2.0
    # Place above fineprint area
    y = fmt.logo_margin
    if fmt.show_fineprint and fmt.fineprint:
        y += fmt.small_fontsize
    parts.append("BT")
    parts.append(f"/F2 {fmt.fontsize} Tf")
    parts.append(f"{x:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(fmt.xaxis_label)}) Tj")
    parts.append("ET")


def _draw_fineprint(parts: list, fmt: LogoFormat, logo_width: int) -> None:
    """Draw fineprint text at bottom-right."""
    assert fmt.line_margin_right is not None
    fp_width = _string_width(fmt.fineprint, fmt.small_fontsize)
    x = logo_width - fmt.logo_margin - fmt.line_margin_right - fp_width
    y = fmt.logo_margin
    parts.append("BT")
    parts.append(f"/F2 {fmt.small_fontsize} Tf")
    parts.append(f"{x:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(fmt.fineprint)}) Tj")
    parts.append("ET")


def _draw_yaxis(parts: list, fmt: LogoFormat, content_x: float,
                content_y: float) -> None:
    """Draw y-axis with tics, labels, and axis label."""
    assert fmt.yaxis_scale is not None

    axis_x = content_x - fmt.stack_margin
    axis_bottom = content_y
    axis_top = content_y + fmt.stack_height

    # Vertical bar
    parts.append("q")
    parts.append(f"{fmt.stroke_width} w")
    parts.append(f"{axis_x:.4f} {axis_bottom:.4f} m")
    parts.append(f"{axis_x:.4f} {axis_top:.4f} l S")
    parts.append("Q")

    # Major tics and labels
    tic_val = 0.0
    while tic_val <= fmt.yaxis_scale + 1e-10:
        y_pos = axis_bottom + (tic_val / fmt.yaxis_scale) * fmt.stack_height

        # Tic mark
        parts.append("q")
        parts.append(f"{fmt.stroke_width} w")
        parts.append(f"{axis_x - fmt.tic_length:.4f} {y_pos:.4f} m")
        parts.append(f"{axis_x:.4f} {y_pos:.4f} l S")
        parts.append("Q")

        # Tic label
        label_text = (
            str(int(tic_val)) if tic_val == int(tic_val) else f"{tic_val:.1f}"
        )
        label_width = _string_width(label_text, fmt.number_fontsize)
        lx = axis_x - fmt.tic_length - fmt.stack_margin - label_width
        ly = y_pos - fmt.number_fontsize * 0.35
        parts.append("BT")
        parts.append(f"/F2 {fmt.number_fontsize} Tf")
        parts.append(f"{lx:.4f} {ly:.4f} Td")
        parts.append(f"({_pdf_escape(label_text)}) Tj")
        parts.append("ET")

        tic_val += fmt.yaxis_tic_interval

    # Minor tics
    assert fmt.yaxis_minor_tic_interval is not None
    if fmt.yaxis_minor_tic_interval > 0:
        minor_val = 0.0
        while minor_val <= fmt.yaxis_scale + 1e-10:
            y_pos = axis_bottom + (minor_val / fmt.yaxis_scale) * fmt.stack_height
            parts.append("q")
            parts.append(f"{fmt.stroke_width} w")
            parts.append(f"{axis_x - fmt.tic_length / 2:.4f} {y_pos:.4f} m")
            parts.append(f"{axis_x:.4f} {y_pos:.4f} l S")
            parts.append("Q")
            minor_val += fmt.yaxis_minor_tic_interval

    # Y-axis label (rotated)
    if fmt.show_yaxis_label and fmt.yaxis_label:
        max_label_str = (
            str(int(fmt.yaxis_scale))
            if fmt.yaxis_scale == int(fmt.yaxis_scale)
            else f"{fmt.yaxis_scale:.1f}"
        )
        label_offset = (
            len(max_label_str) * fmt.number_fontsize * 0.6
            + fmt.tic_length * 1.25
        )
        label_x = axis_x - label_offset
        label_y = axis_bottom + fmt.stack_height / 2

        label_width = _string_width(fmt.yaxis_label, fmt.fontsize)
        # Rotate 90° CCW: the text baseline becomes vertical.
        # We translate to the label position, rotate, then offset by half the
        # text width to center it.
        parts.append("BT")
        parts.append(f"/F2 {fmt.fontsize} Tf")
        # PDF text matrix: [cos sin -sin cos tx ty] for rotation
        # 90° CCW: cos=0, sin=1 → matrix [0 1 -1 0 tx ty]
        parts.append(
            f"0 1 -1 0 {label_x:.4f} {label_y - label_width / 2:.4f} Tm"
        )
        parts.append(f"({_pdf_escape(fmt.yaxis_label)}) Tj")
        parts.append("ET")


def _draw_xaxis_tic(parts: list, fmt: LogoFormat, stack_x: float,
                    stack_y: float, annotation: str) -> None:
    """Draw x-axis tic mark and number annotation below a stack."""
    # Horizontal line at bottom of stack
    parts.append("q")
    parts.append(f"{fmt.stroke_width} w")
    parts.append(f"{stack_x:.4f} {stack_y:.4f} m")
    parts.append(f"{stack_x + fmt.stack_width:.4f} {stack_y:.4f} l S")
    parts.append("Q")

    # Vertical tic
    if annotation:
        tic_len = fmt.tic_length / 2
    else:
        tic_len = fmt.tic_length / 4

    center_x = stack_x + fmt.stack_width / 2
    parts.append("q")
    parts.append(f"{fmt.stroke_width} w")
    parts.append(f"{center_x:.4f} {stack_y:.4f} m")
    parts.append(f"{center_x:.4f} {stack_y - tic_len:.4f} l S")
    parts.append("Q")

    # Number annotation
    if annotation:
        if fmt.rotate_numbers:
            # Rotated 90° CW in PDF coords (text reads top-to-bottom)
            tx = center_x + fmt.number_fontsize * 0.35
            ty = stack_y - tic_len - fmt.stack_margin
            parts.append("BT")
            parts.append(f"/F2 {fmt.number_fontsize} Tf")
            # 90° CW rotation: [0 -1 1 0 tx ty]
            ann_width = _string_width(annotation, fmt.number_fontsize)
            parts.append(
                f"0 -1 1 0 {tx:.4f} {ty:.4f} Tm"
            )
            parts.append(f"({_pdf_escape(annotation)}) Tj")
            parts.append("ET")
        else:
            ann_width = _string_width(annotation, fmt.number_fontsize)
            ax = center_x - ann_width / 2
            ay = stack_y - tic_len - fmt.number_fontsize
            parts.append("BT")
            parts.append(f"/F2 {fmt.number_fontsize} Tf")
            parts.append(f"{ax:.4f} {ay:.4f} Td")
            parts.append(f"({_pdf_escape(annotation)}) Tj")
            parts.append("ET")


def _draw_left_end(parts: list, fmt: LogoFormat, content_x: float,
                   content_y: float) -> None:
    """Draw 5'/N label at left end."""
    if fmt.end_type == "d":
        label = "5"
        use_prime = True
    elif fmt.end_type == "p":
        label = "N"
        use_prime = False
    else:
        return  # pragma: no cover — only called when show_ends is True

    x = content_x - fmt.fontsize * 0.25
    y = content_y - fmt.fontsize * 1.25
    label_width = _string_width(label, fmt.fontsize)

    parts.append("BT")
    parts.append(f"/F2 {fmt.fontsize} Tf")
    parts.append(f"{x - label_width:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(label)}) Tj")
    parts.append("ET")

    if use_prime:
        # Draw prime symbol using Symbol font
        parts.append("BT")
        parts.append(f"/F4 {fmt.fontsize} Tf")
        parts.append(f"{x - label_width + _string_width(label, fmt.fontsize):.4f} {y:.4f} Td")
        parts.append("(\\242) Tj")  # prime in Symbol font
        parts.append("ET")


def _draw_right_end(parts: list, fmt: LogoFormat, right_x: float,
                    content_y: float) -> None:
    """Draw 3'/C label at right end."""
    if fmt.end_type == "d":
        label = "3"
        use_prime = True
    elif fmt.end_type == "p":
        label = "C"
        use_prime = False
    else:
        return  # pragma: no cover — only called when show_ends is True

    x = right_x + fmt.fontsize * 0.25
    y = content_y - fmt.fontsize * 1.25

    parts.append("BT")
    parts.append(f"/F2 {fmt.fontsize} Tf")
    parts.append(f"{x:.4f} {y:.4f} Td")
    parts.append(f"({_pdf_escape(label)}) Tj")
    parts.append("ET")

    if use_prime:
        parts.append("BT")
        parts.append(f"/F4 {fmt.fontsize} Tf")
        parts.append(f"{x + _string_width(label, fmt.fontsize):.4f} {y:.4f} Td")
        parts.append("(\\242) Tj")
        parts.append("ET")


def _draw_symbol(parts: list, symbol: str, x: float, y_bottom: float,
                 target_width: float, target_height: float,
                 color: Color) -> None:
    """Draw a single character scaled to fill the target box.

    In PDF coords, y_bottom is the bottom of the target box, y increases upward.
    """
    if target_height <= 0 or target_width <= 0:
        return  # pragma: no cover — defensive guard

    # IReplacementHack: Helvetica 'I' is too narrow (no serifs).
    # Render two 'T' glyphs overlaid — one normal (crossbar at top) and one
    # flipped vertically (crossbar at bottom) — to create a serifed I shape.
    if symbol == "I":
        _draw_serifed_I(parts, x, y_bottom, target_width, target_height, color)
        return

    ref_size = 100.0
    bbox = _CHAR_BBOX.get(symbol.upper(), _DEFAULT_BBOX)

    lx, ly, ux, uy = bbox
    ref_char_width = (ux - lx) * ref_size
    ref_char_height = (uy - ly) * ref_size

    if ref_char_height <= 0 or ref_char_width <= 0:
        return  # pragma: no cover — all known bboxes have positive dimensions

    sx = target_width / ref_char_width
    sy = target_height / ref_char_height

    tx = x - sx * lx * ref_size
    ty = y_bottom - sy * ly * ref_size

    parts.append("q")
    parts.append(f"{color.red:.4f} {color.green:.4f} {color.blue:.4f} rg")
    parts.append(f"{sx:.6f} 0 0 {sy:.6f} {tx:.4f} {ty:.4f} cm")
    parts.append("BT")
    parts.append(f"/F1 {ref_size} Tf")
    parts.append("0 0 Td")
    parts.append(f"({_pdf_escape(symbol)}) Tj")
    parts.append("ET")
    parts.append("Q")


def _draw_serifed_I(parts: list, x: float, y_bottom: float,
                    target_width: float, target_height: float,
                    color: Color) -> None:
    """Draw a serifed I by overlaying two T glyphs — one upright, one flipped."""
    # Inset slightly to prevent the bottom crossbar from crowding the letter below
    pad = target_height * 0.03
    y_bottom = y_bottom + pad
    target_height = target_height - 2 * pad

    ref_size = 100.0
    bbox = _CHAR_BBOX["T"]
    lx, ly, ux, uy = bbox
    ref_char_width = (ux - lx) * ref_size
    ref_char_height = (uy - ly) * ref_size

    if ref_char_height <= 0 or ref_char_width <= 0:
        return  # pragma: no cover — T bbox always has positive dimensions

    sx = target_width / ref_char_width
    sy = target_height / ref_char_height

    # Normal T — crossbar at top
    tx = x - sx * lx * ref_size
    ty = y_bottom - sy * ly * ref_size

    parts.append("q")
    parts.append(f"{color.red:.4f} {color.green:.4f} {color.blue:.4f} rg")
    parts.append(f"{sx:.6f} 0 0 {sy:.6f} {tx:.4f} {ty:.4f} cm")
    parts.append("BT")
    parts.append(f"/F1 {ref_size} Tf")
    parts.append("0 0 Td")
    parts.append("(T) Tj")
    parts.append("ET")
    parts.append("Q")

    # Flipped T — crossbar at bottom (negative sy, translate to top of box)
    ty_flip = y_bottom + target_height + sy * ly * ref_size

    parts.append("q")
    parts.append(f"{color.red:.4f} {color.green:.4f} {color.blue:.4f} rg")
    parts.append(f"{sx:.6f} 0 0 {-sy:.6f} {tx:.4f} {ty_flip:.4f} cm")
    parts.append("BT")
    parts.append(f"/F1 {ref_size} Tf")
    parts.append("0 0 Td")
    parts.append("(T) Tj")
    parts.append("ET")
    parts.append("Q")


def _draw_errorbar(parts: list, fmt: LogoFormat, stack_x: float,
                   stack_y: float, center: float, down: float,
                   up: float) -> None:
    """Draw error bar on top of a stack."""
    if not fmt.show_errorbars:
        return

    assert fmt.yaxis_scale is not None
    assert fmt.char_width is not None
    points_per_unit = fmt.stack_height / fmt.yaxis_scale
    center_pdf_y = stack_y + center * points_per_unit
    down_y = center_pdf_y - down * points_per_unit * fmt.errorbar_fraction
    up_y = center_pdf_y + up * points_per_unit * fmt.errorbar_fraction

    bar_x = stack_x + fmt.stack_width / 2
    errorbar_width = fmt.char_width * fmt.errorbar_width_fraction
    half_w = errorbar_width / 2

    gray = fmt.errorbar_gray

    parts.append("q")
    parts.append(f"{gray:.4f} {gray:.4f} {gray:.4f} RG")
    parts.append(f"{fmt.stroke_width} w")

    # Vertical line
    parts.append(f"{bar_x:.4f} {down_y:.4f} m")
    parts.append(f"{bar_x:.4f} {up_y:.4f} l S")

    # Bottom cap
    parts.append(f"{bar_x - half_w:.4f} {down_y:.4f} m")
    parts.append(f"{bar_x + half_w:.4f} {down_y:.4f} l S")

    # Top cap
    parts.append(f"{bar_x - half_w:.4f} {up_y:.4f} m")
    parts.append(f"{bar_x + half_w:.4f} {up_y:.4f} l S")

    parts.append("Q")


# ---------------------------------------------------------------------------
# PDF assembly
# ---------------------------------------------------------------------------

def _assemble_pdf(content_stream: str, width: int, height: int) -> bytes:
    """Build a complete, valid PDF 1.4 document."""
    # We build the PDF with 8 objects:
    # 1: Catalog, 2: Pages, 3: Page, 4: Stream (content),
    # 5: Font Helvetica-Bold (F1), 6: Font Helvetica (F2),
    # 7: Font Courier (F3), 8: Font Symbol (F4)

    stream_bytes = content_stream.encode("latin-1")
    stream_length = len(stream_bytes)

    objects: list = []
    offsets: list = []

    header = b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"

    def add_obj(obj_bytes: bytes) -> None:
        offsets.append(len(header) + sum(len(o) for o in objects))
        objects.append(obj_bytes)

    # Object 1: Catalog
    add_obj(b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n")

    # Object 2: Pages
    add_obj(b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n")

    # Object 3: Page
    page = (
        f"3 0 obj\n"
        f"<< /Type /Page /Parent 2 0 R\n"
        f"   /MediaBox [0 0 {width} {height}]\n"
        f"   /Contents 4 0 R\n"
        f"   /Resources <<\n"
        f"     /Font <<\n"
        f"       /F1 5 0 R /F2 6 0 R /F3 7 0 R /F4 8 0 R\n"
        f"     >>\n"
        f"   >>\n"
        f">>\n"
        f"endobj\n"
    )
    add_obj(page.encode("latin-1"))

    # Object 4: Content stream
    stream_obj = (
        f"4 0 obj\n"
        f"<< /Length {stream_length} >>\n"
        f"stream\n"
    ).encode("latin-1") + stream_bytes + b"\nendstream\nendobj\n"
    add_obj(stream_obj)

    # Object 5: Font Helvetica-Bold (F1 — logo characters)
    add_obj(
        b"5 0 obj\n"
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>\n"
        b"endobj\n"
    )

    # Object 6: Font Helvetica (F2 — text labels, numbers)
    add_obj(
        b"6 0 obj\n"
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\n"
        b"endobj\n"
    )

    # Object 7: Font Times-Roman (F3 — I replacement hack, serifed I)
    add_obj(
        b"7 0 obj\n"
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Times-Roman >>\n"
        b"endobj\n"
    )

    # Object 8: Font Symbol (F4 — prime symbol)
    add_obj(
        b"8 0 obj\n"
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Symbol >>\n"
        b"endobj\n"
    )

    # Build body
    body = header + b"".join(objects)

    # Cross-reference table
    xref_offset = len(body)
    num_objects = len(objects) + 1  # +1 for the free entry (object 0)
    xref = f"xref\n0 {num_objects}\n"
    xref += "0000000000 65535 f \n"
    for off in offsets:
        xref += f"{off:010d} 00000 n \n"

    # Trailer
    trailer = (
        f"trailer\n"
        f"<< /Size {num_objects} /Root 1 0 R >>\n"
        f"startxref\n"
        f"{xref_offset}\n"
        f"%%EOF\n"
    )

    return body + xref.encode("latin-1") + trailer.encode("latin-1")
