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

# Funzione per calcolare pI
def calcola_pI(seq):
    try:
        prot = ProteinAnalysis(seq)
        ip = IsoelectricPoint(seq)
        return round(ip.pi(), 2)
    except Exception:
        return "Errore nel calcolo"

# Funzione per contare aromatici
def conta_aromatici(seq):
    aromatici = ['F', 'W', 'Y']
    conteggio = Counter(seq)
    return {aa: conteggio.get(aa, 0) for aa in aromatici}

# Evidenzia aromatici
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
    return f'<a href="data:file/txt;base64,{b64}" download="{filename}">📥 Scarica Sequenza</a>'

# Visualizzazione 3D
def mostra_3D(pdb_id):
    xyzview = py3Dmol.view(query=f"pdb:{pdb_id}")
    xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview, height=400, width=700)

# === Layout ===
col1, col2 = st.columns([1, 1])

with col1:
    st.header("🔎 Cerca proteina")
    nome_proteina = st.text_input("Nome proteina (es. P69905, HBB_HUMAN)", "")
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
            st.success(f"✅ Proteina trovata: {record_id}")
        else:
            st.error("❌ Proteina non trovata. Controlla l'ID UniProt.")

with col2:
    st.header("📁 Oppure carica un file FASTA")
    file = st.file_uploader("Carica file FASTA", type=["fasta", "fa", "txt"])
    if file:
        try:
            contenuto = file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(contenuto), "fasta")
            sequenza = str(record.seq)
            record_id = record.id
            st.success(f"✅ File caricato: {record_id}")
        except Exception:
            st.error("❌ Errore nel caricamento del file")

# === Analisi ===
if sequenza:
    with st.expander("📊 Analisi della sequenza", expanded=True):
        st.markdown(f"**ID:** `{record_id}`")
        st.markdown(f"**Lunghezza:** {len(sequenza)} amminoacidi")
        st.markdown("**Sequenza:**")
        st.markdown(evidenzia_aromatici(sequenza), unsafe_allow_html=True)

        st.markdown("**Aromatici (F, W, Y):**")
        st.write(conta_aromatici(sequenza))

        st.markdown(f"**Punto isoelettrico (pI):** {calcola_pI(sequenza)}")

        st.markdown(download_link(sequenza, "sequenza_proteica.txt"), unsafe_allow_html=True)

        # Grafico
        st.markdown("**🧪 Frequenza amminoacidi:**")
        conta = Counter(sequenza)
        fig, ax = plt.subplots()
        ax.bar(conta.keys(), conta.values(), color="skyblue")
        plt.xlabel("Amminoacidi")
        plt.ylabel("Frequenza")
        st.pyplot(fig)

    # === Visualizzazione 3D ===
    with st.expander("🧬 Visualizzazione 3D (solo per proteine note con PDB)", expanded=False):
        pdb_id = st.text_input("Inserisci PDB ID (es. 1A3N, 4HHB)")
        if pdb_id:
            try:
                mostra_3D(pdb_id)
            except Exception as e:
                st.error(f"❌ Errore visualizzazione: {e}")
else:
    st.info("Inserisci un nome proteina oppure carica un file FASTA per iniziare.")

st.markdown(
    """
    <div style="text-align: center; margin-top: 50px;">
        <a href="https://www.buymeacoffee.com/gjovannj" target="_blank">
            <img src="https://img.buymeacoffee.com/button-api/?text=Offrimi un caffè&emoji=☕&slug=gjovannj&button_colour=FFDD00&font_colour=000000&font_family=Arial&outline_colour=000000&coffee_colour=ffffff" />
        </a>
    </div>
    """,
    unsafe_allow_html=True
)



