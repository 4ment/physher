"""The Streamlit front end: ``phyforge gui``.

Like the CLI it only turns input into an :class:`~phyforge.spec.AnalysisSpec`;
the builders do the rest.  It needs the ``gui`` extra (streamlit, plotly,
pandas, numpy), which nothing outside this package imports.
"""
