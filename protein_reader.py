import streamlit as st
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import io
import base64

st.set_page_config(page_title="Protein Explorer", layout="wide")

def evidenzia_aromatici(seq):
    # Evidenzia F, W, Y con colore rosso
    html_seq = ""
    for aa in seq:
        if aa in ("F", "W", "Y"):
            html_seq += f'<span style="color: red; font-weight: bold;">{aa}</span>'
        else:
            html_seq += aa
    return html_seq

def calcola_pI(seq):
    # Calcolo semplificato punto isoelettrico (valore indicativo)
    # Conta amminoacidi basici (K,R,H) e acidi (D,E)
    basici = seq.count('K') + seq.count('R') + seq.count('H')
    acidi = seq.count('D') + seq.count('E')
    # pI approssimato: più basici -> pI alto, più acidi -> pI basso
    return 7 + (basici - acidi)*0.1

def crea_link_download(text, filename, label):
    b64 = base64.b64encode(text.encode()).decode()  # convert to base64
    href = f'<a href="data:file/txt;base64,{b64}" download="{filename}">{label}</a>'
    return href

st.title("Protein Explorer")

file_caricato = st.file_uploader("Carica un file FASTA (supporta più sequenze)", type=["fasta", "fa", "txt"])

if file_caricato is not None:
    try:
        contenuto_testo = io.StringIO(file_caricato.getvalue().decode("utf-8"))
        records = list(SeqIO.parse(contenuto_testo, "fasta"))
        if len(records) == 0:
            st.error("Nessuna sequenza trovata nel file FASTA.")
        else:
            st.success(f"Trovate {len(records)} sequenze nel file.")

            for i, record in enumerate(records):
                with st.expander(f"Sequenza {i+1}: {record.id}"):
                    seq = str(record.seq)

                    col1, col2 = st.columns([3,2])

                    with col1:
                        st.subheader("Sequenza proteica:")
                        st.markdown(evidenzia_aromatici(seq), unsafe_allow_html=True)

                    with col2:
                        st.subheader("Analisi")
                        counts = Counter(seq)
                        aromatici = {aa: counts.get(aa, 0) for aa in ("F", "W", "Y")}
                        st.write("Conteggio amminoacidi aromatici (F, W, Y):")
                        st.write(aromatici)

                        pI = calcola_pI(seq)
                        st.write(f"Punto isoelettrico (pI) approssimato: **{pI:.2f}**")

                        # Grafico frequenze
                        st.write("Frequenza amminoacidi:")
                        fig, ax = plt.subplots(figsize=(6,3))
                        aminoacidi = list(counts.keys())
                        frequenze = list(counts.values())
                        ax.bar(aminoacidi, frequenze, color='skyblue')
                        ax.set_xlabel("Amminoacidi")
                        ax.set_ylabel("Frequenza")
                        ax.set_title("Conteggio amminoacidi")
                        st.pyplot(fig)

                    # Link per download sequenza e analisi
                    testo_download = f"> {record.description}\n{seq}\n\nConteggio aromatici: {aromatici}\nPunto isoelettrico: {pI:.2f}"
                    st.markdown(crea_link_download(testo_download, f"analisi_{record.id}.txt", "📥 Scarica analisi"), unsafe_allow_html=True)

    except Exception as e:
        st.error(f"Errore nel caricamento o parsing del file: {e}")

else:
    st.info("Carica un file FASTA per iniziare l'analisi.")


