html = """
<!DOCTYPE html>
<html>
    <head>
        <meta charset='utf-8'>
        <title>NASA GeneLab Data API reference</title>
        <style type='text/css'>
            span.code {{
                border: 1px solid #888888;
                padding: 5pt;
            }}
            span.code a, span.code tt {{
                font-family: monospace, mono;
            }}
        </style>
    </head>

    <body>
        <h1>NASA GeneLab Data API reference</h1>
            <h2>The anatomy of a request</h2>
                <span class='code'>
                    <tt>{url_root}</tt> <tt>/</tt>
                    <a style='background:#eedd55;'>dataset_accession</a> <tt>/</tt>
                    <a style='background:#55dd55;'>assay_name</a> <tt>/</tt>
                    <a style='background:#dd5555;'>data_category</a> <tt>/</tt>
                    <a style='background:#aaaaff;'>data_type</a> <tt>/</tt>
                    <a style='background:#ff88ff;'>transform</a> <tt>/?</tt>
                    <a style='background:#aaaaaa;'>[get_args]</a>
                </span>
    </body>
</html>
"""
