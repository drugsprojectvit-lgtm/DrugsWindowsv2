import base64
import json
from IPython.display import HTML

def show_structure(protein_text: str = None, ligand_text: str = None, pdb_id: str = "Structure", protein_name: str = "") -> str:
    """
    Robust 3D visualization supporting:
    1. Protein Only
    2. Protein + Ligand (Overlay)
    3. Ligand Only (Standalone)
    """
    
    # 1. Safely serialize inputs
    # If None is passed, json.dumps(None) produces the string "null" (without quotes in JS), which is perfect.
    prot_json = json.dumps(protein_text)
    lig_json = json.dumps(ligand_text)
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {{ margin: 0; padding: 0; overflow: hidden; background-color: white; }}
            #container {{ width: 100vw; height: 100vh; position: relative; }}
            #error-log {{ 
                display: none; position: absolute; top: 10px; left: 10px; 
                background: rgba(255,0,0,0.8); color: white; padding: 10px; 
                z-index: 999; font-family: sans-serif; border-radius: 5px;
            }}
            .legend {{
                position: absolute; bottom: 10px; right: 10px;
                background: rgba(255, 255, 255, 0.9); padding: 8px 12px;
                border-radius: 6px; font-family: sans-serif; font-size: 12px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2); z-index: 100;
                pointer-events: none;
            }}
            .color-box {{ 
                display: inline-block; width: 10px; height: 10px; 
                margin-right: 6px; border-radius: 50%; 
            }}
            .info-overlay {{
                position: absolute; top: 10px; left: 10px;
                background: rgba(255, 255, 255, 0.9); padding: 8px 12px;
                border-radius: 6px; font-family: sans-serif; font-size: 14px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.2); z-index: 90;
            }}
        </style>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="error-log"></div>
        
        <div class="info-overlay">
            <strong>{pdb_id}</strong><br>
            <span style="font-size:12px; color:#555">{protein_name}</span>
        </div>

        <div class="legend">
            <div><span class="color-box" style="background: linear-gradient(90deg, blue, green, red);"></span>Protein</div>
            <div style="margin-top:4px"><span class="color-box" style="background: magenta;"></span>Ligand</div>
        </div>

        <div id="container"></div>

        <script>
            function logError(msg) {{
                var el = document.getElementById('error-log');
                el.style.display = 'block';
                el.innerHTML += "Error: " + msg + "<br>";
                console.error(msg);
            }}

            window.onload = function() {{
                try {{
                    // 1. Check Library
                    if (typeof $3Dmol === 'undefined') {{
                        throw new Error("3Dmol.js failed to load. Check internet connection.");
                    }}

                    // 2. Initialize Viewer
                    var element = document.getElementById('container');
                    var config = {{ backgroundColor: 'white' }};
                    var viewer = $3Dmol.createViewer(element, config);

                    var proteinData = {prot_json};
                    var ligandData = {lig_json};
                    var hasModel = false;

                    // --- ADD PROTEIN (If exists) ---
                    if (proteinData) {{
                        viewer.addModel(proteinData, "pdb");
                        // Style last added model (Protein)
                        viewer.setStyle(
                            {{model: -1}}, 
                            {{cartoon: {{color: 'spectrum'}}}}
                        );
                        // Optional: Show existing hetatoms in protein file as faint sticks
                        viewer.addStyle(
                            {{model: -1, hetflag: true}}, 
                            {{stick: {{radius: 0.1, color: 'lightgray'}}}}
                        );
                        hasModel = true;
                    }}

                    // --- ADD LIGAND (If exists) ---
                    if (ligandData) {{
                        viewer.addModel(ligandData, "pdb");
                        // Style last added model (Ligand)
                        viewer.setStyle(
                            {{model: -1}}, 
                            {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}}
                        );
                        hasModel = true;
                    }} else if (proteinData) {{
                        // Fallback: If no separate ligand file, look for ligands INSIDE protein
                        viewer.addStyle(
                            {{resn: ["UNL", "LIG", "DRG", "UNK"]}}, 
                            {{stick: {{colorscheme: 'magentaCarbon', radius: 0.4}}}}
                        );
                    }}

                    if (!hasModel) {{
                        throw new Error("No protein or ligand data provided.");
                    }}

                    // 4. Render
                    viewer.zoomTo();
                    viewer.render();

                }} catch(e) {{
                    logError(e.message);
                }}
            }};
        </script>
    </body>
    </html>
    """
    
    b64 = base64.b64encode(html_content.encode()).decode()
    iframe = f'<iframe src="data:text/html;base64,{b64}" width="100%" height="500" frameborder="0" style="border: 1px solid #ccc; border-radius: 8px;"></iframe>'
    
    return iframe