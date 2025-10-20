import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="DNA Sequence Feature & Structure Visualizer",
    layout="wide",
    page_icon="🧬"
)

st.title("🧬 DNA Sequence Feature & Structure Visualizer")
st.markdown("""
Paste your nucleotide sequence below to calculate features, k-mers, codons, predict hairpins, and visualize GC content.
""")

# -----------------------------
# Input Sequence
# -----------------------------
seq_input = st.text_area("Enter DNA Sequence (A/T/G/C)", height=200)

if seq_input:
    seq_str = seq_input.replace("\n","").upper()
    length = len(seq_str)

    # -----------------------------
    # Basic Features
    # -----------------------------
    gc_content = 100 * (seq_str.count('G') + seq_str.count('C')) / length
    at_content = 100 * (seq_str.count('A') + seq_str.count('T')) / length
    at_gc_skew = ((seq_str.count('A') + seq_str.count('T')) - (seq_str.count('G') + seq_str.count('C'))) / length

    # -----------------------------
    # Dinucleotide Counts
    # -----------------------------
    dinucleotides = [a+b for a in "ATGC" for b in "ATGC"]
    di_counts = {k: seq_str.count(k) for k in dinucleotides}

    # -----------------------------
    # Trinucleotide Counts
    # -----------------------------
    trinucleotides = [a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"]
    tri_counts = {k: seq_str.count(k) for k in trinucleotides if seq_str.count(k)>0}

    # -----------------------------
    # Codon Usage
    # -----------------------------
    codons = {}
    if length % 3 == 0:
        for i in range(0, length, 3):
            codon = seq_str[i:i+3]
            codons[codon] = codons.get(codon,0)+1

    # -----------------------------
    # Hairpin / Palindrome Finder
    # -----------------------------
    def find_palindromes(seq, min_len=4, max_len=12):
        pals = []
        for l in range(min_len, max_len+1):
            for i in range(len(seq)-l+1):
                fragment = seq[i:i+l]
                revcomp = fragment[::-1].translate(str.maketrans("ATGC", "TACG"))
                if fragment == revcomp:
                    pals.append((fragment, i+1, i+l))
        return pals

    hairpins = find_palindromes(seq_str)

    # -----------------------------
    # GC Sliding Window
    # -----------------------------
    def gc_sliding_window(seq, window=10):
        gc_vals = []
        for i in range(len(seq)-window+1):
            fragment = seq[i:i+window]
            gc = 100*(fragment.count('G')+fragment.count('C'))/window
            gc_vals.append(gc)
        return gc_vals

    gc_values = gc_sliding_window(seq_str, window=10)

    # -----------------------------
    # Folding Sections
    # -----------------------------
    with st.expander("Basic Sequence Features", expanded=True):
        features = {
            "Length": length,
            "GC Content (%)": round(gc_content,2),
            "AT Content (%)": round(at_content,2),
            "AT-GC Skew": round(at_gc_skew,3)
        }
        st.dataframe(pd.DataFrame.from_dict(features, orient='index', columns=['Value']))

    with st.expander("Dinucleotide Counts", expanded=False):
        di_df = pd.DataFrame.from_dict(di_counts, orient='index', columns=['Count'])
        st.dataframe(di_df)
        fig, ax = plt.subplots(figsize=(12,4))
        ax.bar(di_df.index, di_df['Count'], color='skyblue')
        ax.set_ylabel("Count")
        ax.set_title("Dinucleotide Frequencies")
        plt.xticks(rotation=45)
        st.pyplot(fig)

    if tri_counts:
        with st.expander("Trinucleotide Counts", expanded=False):
            tri_df = pd.DataFrame.from_dict(tri_counts, orient='index', columns=['Count'])
            st.dataframe(tri_df)
            fig2, ax2 = plt.subplots(figsize=(12,4))
            ax2.bar(tri_df.index, tri_df['Count'], color='lightgreen')
            ax2.set_ylabel("Count")
            ax2.set_title("Trinucleotide Frequencies")
            plt.xticks(rotation=90)
            st.pyplot(fig2)

    if codons:
        with st.expander("Codon Usage", expanded=False):
            codon_df = pd.DataFrame.from_dict(codons, orient='index', columns=['Count'])
            st.dataframe(codon_df)
            fig3, ax3 = plt.subplots(figsize=(12,4))
            ax3.bar(codon_df.index, codon_df['Count'], color='salmon')
            ax3.set_ylabel("Count")
            ax3.set_title("Codon Usage")
            plt.xticks(rotation=90)
            st.pyplot(fig3)

    if hairpins:
        with st.expander("Potential Hairpin / Palindromic Regions", expanded=False):
            hp_df = pd.DataFrame(hairpins, columns=["Sequence","Start","End"])
            st.dataframe(hp_df)
            st.markdown(f"Found **{len(hairpins)}** potential hairpin motifs.")
    else:
        st.info("No palindromic hairpin motifs detected.")

    with st.expander("GC Content Sliding Window (window=10)", expanded=False):
        fig4, ax4 = plt.subplots(figsize=(12,4))
        ax4.plot(range(len(gc_values)), gc_values, color='purple')
        ax4.set_xlabel("Position")
        ax4.set_ylabel("GC %")
        ax4.set_title("GC Content Across Sequence")
        st.pyplot(fig4)