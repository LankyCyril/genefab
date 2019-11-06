html = """
<!DOCTYPE html>
<html>
    <head>
        <meta charset='utf-8'>
        <title>NASA GeneLab Data API reference</title>
        <style type='text/css'>
            span.code {{
                font-family: monospace, mono;
                border: 1px solid #888888;
                padding: 5pt;
            }}
            span.code a, span.code tt {{
                font-family: monospace, mono;
            }}
            code {{
                font-family: monospace, mono;
                font-size: .85em;
            }}
            code a {{
                color: #000;
                text-decoration: none;
            }}
        </style>
    </head>

    <body>
        <h1>NASA GeneLab Data API reference</h1>
            <h2><a name='request_anatomy'>The anatomy of a request</a></h2>
                <span class='code'>
                    <tt>{url_root}</tt> <tt>/</tt>
                    <a style='background:#eedd55;'>dataset_accession</a> /
                    <a style='background:#55dd55;'>assay_name</a> /
                    <a style='background:#dd5555;'>data_category</a> /
                    <a style='background:#aaaaff;'>data_type</a> /
                    <a style='background:#ff88ff;'>transform</a> /?
                    <a style='background:#cccccc;'>[get_args]</a>
                </span>
                <p>
                    The request can be routed to any depth. All following URLs
                    are valid examples:
                </p>
                <code>
                    <a href='{url_root}/GLDS-111/'>{url_root}/<span style='background:#eedd55'>GLDS-111</span>/</a><br>
                    <a href='{url_root}/GLDS-4/assay/?fmt=html'>{url_root}/<span style='background:#eedd55'>GLDS-4</span>/<span style='background:#55dd55'>assay</span>/?<span style='background:#cccccc'>fmt=html</span></a><br>
                    <a href='{url_root}/GLDS-42/assay/annotation/?fmt=json'>{url_root}/<span style='background:#eedd55'>GLDS-42</span>/<span style='background:#55dd55'>assay</span>/<span style='background:#dd5555'>annotation</span>/?<span style='background:#cccccc'>fmt=json</span></a><br>
                    <a href='{url_root}/GLDS-30/a_GLDS-30_microarray_metadata-txt/data/pca/'>{url_root}/<span style='background:#eedd55'>GLDS-30</span>/<span style='background:#55dd55'>a_GLDS-30_microarray_metadata-txt</span>/<span style='background:#dd5555'>data</span>/<span style='background:#aaaaff'>pca</span>/</a><br>
                    <a href='{url_root}/GLDS-63/assay/data/processed/descriptive/?sort_by=Parameter%20Value:%20Absorbed%20Radiation%20Dose&hidecol=Protocol%20REF&top=20'>{url_root}/<span style='background:#eedd55'>GLDS-63</span>/<span style='background:#55dd55'>assay</span>/<span style='background:#dd5555'>data</span>/<span style='background:#aaaaff'>processed</span>/<span style='background:#ff88ff'>descriptive</span>/?<span style='background:#cccccc'>sort_by=Parameter%20Value:%20Absorbed%20Radiation%20Dose&hidecol=Protocol%20REF&top=20</span></a>
                </code>
    </body>
</html>
"""
