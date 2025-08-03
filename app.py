import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis, IsoelectricPoint
from collections import Counter
import matplotlib.pyplot as plt
from io import StringIO
import requests
from stmol import showmol
import py3Dmol
import base64

st.set_page_config(page_title="Protein Explorer", layout="wide")

st.title("🧬 Protein Explorer")

# Function to calculate pI
def calcola_pI(seq):
    try:
        prot = ProteinAnalysis(seq)
        ip = IsoelectricPoint(seq)
        return round(ip.pi(), 2)
    except Exception:
        return "Error in calculation"

# Function to count aromatics
def conta_aromatici(seq):
    aromatici = ['F', 'W', 'Y']
    conteggio = Counter(seq)
    return {aa: conteggio.get(aa, 0) for aa in aromatici}

# Highlight aromatics
def evidenzia_aromatici(seq):
    evidenziata = ""
    for aa in seq:
        if aa in ["F", "W", "Y"]:
            evidenziata += f"<span style='color:red; font-weight:bold'>{aa}</span>"
        else:
            evidenziata += aa
    return f"<div style='font-family: monospace'>{evidenziata}</div>"

# Download helper
def download_link(text, filename):
    b64 = base64.b64encode(text.encode()).decode()
    return f'<a href="data:file/txt;base64,{b64}" download="{filename}">📥 Download Sequence</a>'

# 3D Visualization
def mostra_3D(pdb_id):
    xyzview = py3Dmol.view(query=f"pdb:{pdb_id}")
    xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview, height=400, width=700)

# === Layout ===
col1, col2 = st.columns([1, 1])

with col1:
    st.header("🔎 Search Protein")
    nome_proteina = st.text_input("Protein name (e.g. P69905, HBB_HUMAN)", "")
    sequenza = ""
    record_id = ""

    if nome_proteina:
        url = f"https://rest.uniprot.org/uniprotkb/{nome_proteina}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            fasta_text = response.text
            record = SeqIO.read(StringIO(fasta_text), "fasta")
            sequenza = str(record.seq)
            record_id = record.id
            st.success(f"✅ Protein found: {record_id}")
        else:
            st.error("❌ Protein not found. Check the UniProt ID.")

with col2:
    st.header("📁 Or Upload FASTA File")
    file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])
    if file:
        try:
            contenuto = file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(contenuto), "fasta")
            sequenza = str(record.seq)
            record_id = record.id
            st.success(f"✅ File uploaded: {record_id}")
        except Exception:
            st.error("❌ Error during file processing")

# === Sequence Analysis ===
if sequenza:
    with st.expander("📊 Sequence Analysis", expanded=True):
        st.markdown(f"**ID:** `{record_id}`")
        st.markdown(f"**Length:** {len(sequenza)} amino acids")
        st.markdown("**Sequence:**")
        st.markdown(evidenzia_aromatici(sequenza), unsafe_allow_html=True)

        st.markdown("**Aromatics (F, W, Y):**")
        st.write(conta_aromatici(sequenza))

        st.markdown(f"**Isoelectric point (pI):** {calcola_pI(sequenza)}")

        st.markdown(download_link(sequenza, "protein_sequence.txt"), unsafe_allow_html=True)

        st.markdown("**🧪 Amino Acid Frequency:**")
        conta = Counter(sequenza)
        fig, ax = plt.subplots()
        ax.bar(conta.keys(), conta.values(), color="skyblue")
        plt.xlabel("Amino Acids")
        plt.ylabel("Frequency")
        st.pyplot(fig)

    with st.expander("🧬 3D Visualization (only for proteins with PDB ID)", expanded=False):
        pdb_id = st.text_input("Enter PDB ID (e.g. 1A3N, 4HHB)")
        if pdb_id:
            try:
                mostra_3D(pdb_id)
            except Exception as e:
                st.error(f"❌ 3D rendering error: {e}")
else:
    st.info("Upload a FASTA file or enter a UniProt ID to begin")

# Buy Me a Coffee button
st.markdown(
    """
    <div style="text-align: center; margin-top: 50px;">
        <a href="https://www.buymeacoffee.com/gjovannj" target="_blank">
            <img src="https://img.buymeacoffee.com/button-api/?text=Buy me a coffee&emoji=☕&slug=gjovannj&button_colour=FFDD00&font_colour=000000&font_family=Arial&outline_colour=000000&coffee_colour=ffffff" />
        </a>
    </div>
    """,
    unsafe_allow_html=True
)
