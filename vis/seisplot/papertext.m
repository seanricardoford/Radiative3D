# papertext.m
#
# Like the text() graphics function, but plots in paperspace rather
# than within an axes. Usage:
#
#   h = papertext(X, Y, "text", [prop, val], ... )
#
# X and Y both range from 0 to 1 and span the entire paperspace of the
# figure.
#
# Uses a figure annotation object so the text is positioned in normalized
# figure coordinates without adding a second axes.  The gnuplot toolkit
# renders multiple axes as a multiplot and emits warnings for its streamed
# data when printing, so keeping annotations on the figure avoids that path.
#
function h = papertext(varargin)

  if (length(varargin) < 3)
     error("Usage: h = papertext(X, Y, \"text\", [prop, val], ... )");
  end

  x = varargin{1};
  y = varargin{2};
  label = varargin{3};
  properties = varargin(4:end);

  h = annotation(gcf(), "textbox", [x, y, 0, 0],
                 "string", label,
                 "fitboxtotext", "on",
                 "linestyle", "none",
                 properties{:});

end
