"""Self-contained HTML/JS embeds for the alignment and tree viewers.

Ported from torchseed so the two front ends look the same.  Each function
returns a complete HTML document for ``streamlit.components.v1.html``; the
libraries come from a CDN, tried in turn.
"""

from __future__ import annotations

import json
from typing import Sequence, Tuple

__all__ = ["msaviewer_html", "phylocanvas_html"]


def msaviewer_html(
    sequences: Sequence[Tuple[str, str]], protein: bool = False, height: int = 600
) -> str:
    """Return an HTML document that renders the alignment in BioJS MSAViewer."""
    fasta_text = "\n".join(f">{name}\n{sequence}" for name, sequence in sequences)
    fasta_js = json.dumps(fasta_text)  # safely escaped string literal for JS
    palette = "taylor" if protein else "nucleotide"
    align_h = max(height - 80, 200)
    # NOTE: the `*.min.gz.js` builds shipped in the npm package are *raw gzip
    # bytes* served without `Content-Encoding: gzip`, so browsers can't execute
    # them. Use the uncompressed bundle instead.
    cdn_urls = [
        "https://cdn.jsdelivr.net/npm/msa@1.0.3/dist/msa.js",
        "https://unpkg.com/msa@1.0.3/dist/msa.js",
    ]
    cdn_js = json.dumps(cdn_urls)
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <style>
    html, body {{ margin: 0; padding: 0; font-family: sans-serif; }}
    #status {{ padding: 6px 10px; font-size: 12px; color: #555; }}
    #err    {{ padding: 6px 10px; color: #b00; font-family: monospace;
              white-space: pre-wrap; }}
    #msa    {{ width: 100%; }}
  </style>
</head>
<body>
  <div id="status">Loading BioJS MSAViewer…</div>
  <div id="err"></div>
  <div id="msa"></div>
  <script>
    (function() {{
      var statusEl = document.getElementById("status");
      var errEl    = document.getElementById("err");
      var urls     = {cdn_js};
      var fastaText = {fasta_js};

      window.addEventListener("error", function(e) {{
        errEl.textContent += "[error] " + (e.message || e) + "\\n";
      }});

      // ES2022 added Array.prototype.at, which breaks msa@1.0.3's
      // `"at" in t` branch-detection in resetSeqs. Strip the offending
      // statement from the bundle source before executing it.
      // https://github.com/wilzbach/msa/issues/257
      var BAD = 't.on("add change reset ",e,this),e.call(this)';

      function inject(text) {{
        text = text.replace(BAD, "");
        var s = document.createElement("script");
        s.text = text;
        document.head.appendChild(s);
      }}

      function tryLoad(i) {{
        if (i >= urls.length) {{
          statusEl.style.display = "none";
          errEl.textContent = "Failed to load BioJS MSAViewer from any CDN. "
            + "Tried:\\n  " + urls.join("\\n  ");
          return;
        }}
        statusEl.textContent = "Fetching " + urls[i] + " …";
        fetch(urls[i])
          .then(function(r) {{ if (!r.ok) throw new Error(r.status); return r.text(); }})
          .then(function(text) {{
            inject(text);
            if (typeof msa === "undefined") {{ tryLoad(i + 1); return; }}
            statusEl.style.display = "none";
            try {{
              var m = new msa({{
                el: document.getElementById("msa"),
                colorscheme: {{ scheme: "{palette}" }},
                vis: {{ conserv: false, overviewbox: false, seqlogo: false,
                       labelId: false, metacell: false }},
                zoomer: {{ alignmentHeight: {align_h}, labelNameLength: 150 }},
              }});
              m.seqs.reset(msa.io.fasta.parse(fastaText));
              m.render();
            }} catch (e) {{
              errEl.textContent = "MSAViewer init error: " + e.message;
            }}
          }})
          .catch(function() {{ tryLoad(i + 1); }});
      }}
      tryLoad(0);
    }})();
  </script>
</body>
</html>"""


def phylocanvas_html(newick: str, height: int = 600) -> str:
    """Return an HTML document that renders a Newick tree with phylocanvas.gl.

    Toolbar exposes zoom in/out, vertical leaf-spacing +/-, layout switch, and reset.
    """
    newick_js = json.dumps(newick)
    cdn_urls = [
        "https://cdn.jsdelivr.net/npm/@phylocanvas/phylocanvas.gl@1.63.0/dist/bundle.min.js",
        "https://unpkg.com/@phylocanvas/phylocanvas.gl@1.63.0/dist/bundle.min.js",
    ]
    cdn_js = json.dumps(cdn_urls)
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <style>
    html, body {{ margin: 0; padding: 0; font-family: sans-serif; height: 100%; }}
    #toolbar {{ display: flex; flex-wrap: wrap; gap: 4px; padding: 6px 8px;
               border-bottom: 1px solid #ddd; background: #f7f7f7;
               font-size: 12px; align-items: center; }}
    #toolbar button {{ padding: 3px 8px; font-size: 12px; cursor: pointer;
                      border: 1px solid #bbb; background: #fff; border-radius: 3px; }}
    #toolbar button:hover {{ background: #eee; }}
    #toolbar select {{ font-size: 12px; padding: 2px 4px; }}
    #toolbar .sep {{ width: 1px; height: 18px; background: #ccc; margin: 0 4px; }}
    #toolbar .label {{ color: #555; }}
    #err {{ padding: 8px; color: #b00; font-family: monospace;
           white-space: pre-wrap; }}
    #tree {{ width: 100%; height: calc(100% - 40px); position: relative; }}
  </style>
</head>
<body>
  <div id="toolbar">
    <button id="zin"  title="Zoom in">＋ Zoom</button>
    <button id="zout" title="Zoom out">－ Zoom</button>
    <span class="sep"></span>
    <span class="label">Spacing:</span>
    <button id="vin"  title="Increase leaf spacing">↕＋</button>
    <button id="vout" title="Decrease leaf spacing">↕−</button>
    <span class="sep"></span>
    <span class="label">Layout:</span>
    <select id="layout">
      <option value="rc">Rectangular</option>
      <option value="dg">Diagonal</option>
      <option value="hr">Hierarchical</option>
      <option value="rd">Radial</option>
      <option value="cr">Circular</option>
    </select>
    <span class="sep"></span>
    <button id="reset" title="Reset view">Reset</button>
  </div>
  <div id="err"></div>
  <div id="tree"></div>
  <script>
    (function() {{
      var errEl   = document.getElementById("err");
      var treeEl  = document.getElementById("tree");
      var urls    = {cdn_js};
      var newick  = {newick_js};

      window.addEventListener("error", function(e) {{
        errEl.textContent += "[error] " + (e.message || e) + "\\n";
      }});

      function tryLoad(i) {{
        if (i >= urls.length) {{
          errEl.textContent = "Failed to load phylocanvas.gl from any CDN. Tried:\\n  "
            + urls.join("\\n  ");
          return;
        }}
        var s = document.createElement("script");
        s.src = urls[i];
        s.onload = function() {{
          if (typeof phylocanvas === "undefined") {{ tryLoad(i + 1); return; }}
          init();
        }};
        s.onerror = function() {{ tryLoad(i + 1); }};
        document.head.appendChild(s);
      }}

      function init() {{
        try {{
          var PG   = phylocanvas.PhylocanvasGL;
          var TT   = phylocanvas.TreeTypes || {{
            Rectangular: "rc", Diagonal: "dg", Hierarchical: "hr",
            Radial: "rd", Circular: "cr",
          }};
          var rect = treeEl.getBoundingClientRect();
          var tree = new PG(treeEl, {{
            size: {{ width: rect.width, height: rect.height }},
            source: newick,
            type: TT.Rectangular,
            showLabels: true,
            nodeSize: 1,
            showLeafLabels: true,
            interactive: true,
            padding: 20,
          }});

          // Defaults to remember for Reset
          var def = {{
            zoom:         tree.getZoom ? tree.getZoom() : 0,
            stepZoom:     0,
            branchScalar: 1,
            type:         TT.Rectangular,
          }};

          function bump(prop, delta, mult) {{
            var cur = tree.props && tree.props[prop];
            if (cur === undefined || cur === null) cur = (prop === "branchScalar" ? 1 : 0);
            var next = mult ? cur * delta : cur + delta;
            var p = {{}}; p[prop] = next;
            tree.setProps(p);
          }}

          document.getElementById("zin").onclick  = function() {{ bump("zoom",  0.5, false); }};
          document.getElementById("zout").onclick = function() {{ bump("zoom", -0.5, false); }};
          document.getElementById("vin").onclick  = function() {{ bump("stepZoom",  0.25, false); }};
          document.getElementById("vout").onclick = function() {{ bump("stepZoom", -0.25, false); }};
          document.getElementById("layout").onchange = function(ev) {{
            tree.setProps({{ type: ev.target.value }});
          }};
          document.getElementById("reset").onclick = function() {{
            tree.setProps({{
              zoom: def.zoom, stepZoom: def.stepZoom,
              branchScalar: def.branchScalar, type: def.type,
            }});
            document.getElementById("layout").value = def.type;
          }};

          window.addEventListener("resize", function() {{
            var r = treeEl.getBoundingClientRect();
            tree.setProps({{ size: {{ width: r.width, height: r.height }} }});
          }});
        }} catch (e) {{
          errEl.textContent = "Phylocanvas init error: " + (e && e.message ? e.message : e);
        }}
      }}
      tryLoad(0);
    }})();
  </script>
</body>
</html>"""
