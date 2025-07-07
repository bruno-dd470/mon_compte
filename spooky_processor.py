# -*- coding: utf-8 -*-
# =============================================================================
# Standard Library Imports
# =============================================================================
import os
import sys
import re
import logging
import time
import csv
import threading
import platform
import socket  # <-- Ajout
import random  # <-- Ajout
import json    # <-- Ajout
import math
from datetime import datetime
from pathlib import Path
from collections import OrderedDict
from urllib.error import HTTPError
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError

# =============================================================================
# Configuration Globale
# =============================================================================
# Sécurité SSL (désactive la vérification pour les requêtes NCBI)
import ssl
ssl._create_default_https_context = ssl._create_unverified_context  # <-- Après import ssl

# Timeout réseau global (60 secondes)
socket.setdefaulttimeout(60)  # <-- Après import socket

# Verrou pour les logs thread-safe
log_lock = threading.Lock()  # <-- Après import threading

# =============================================================================
# Third-Party Imports
# =============================================================================
import numpy as np
import requests
import gzip
from Bio import Entrez, SeqIO
from io import BytesIO, StringIO

# =============================================================================
# GUI Imports (conditionnels)
# =============================================================================
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, simpledialog
except ImportError:
    pass  # Mode CLI seulement
def get_config_path():
    """Retourne le chemin optimal selon l'OS"""
    system = platform.system()
    
    if system == "Windows":
        base_dir = Path(os.environ.get('APPDATA'))
    elif system == "Darwin":  # macOS
        base_dir = Path.home() / "Library" / "Application Support"
    else:  # Linux/Unix
        base_dir = Path.home() / ".config"
    
    app_dir = base_dir / "SpookyProcessor"
    app_dir.mkdir(exist_ok=True)
    return app_dir / "ncbi_config.json"

def load_ncbi_config():
    config_path = get_config_path()
    if config_path.exists():
        try:
            with open(config_path, "r") as f:
                return json.load(f)
        except Exception as e:
            logging.warning(f"Invalid config file: {str(e)}")
    return {}

def save_ncbi_config(email, api_key=None):
    config_path = get_config_path()
    config = {"email": email}
    if api_key:
        config["api_key"] = api_key
    
    try:
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)
        os.chmod(config_path, 0o600)  # Permissions restrictives (Unix)
    except Exception as e:
        logging.error(f"Failed to save config: {str(e)}")

# Configuration NCBI par défaut
def init_ncbi():
    config = load_ncbi_config()
    Entrez.email = config.get("email", "your_email@example.com")  # Email par défaut obligatoire
    Entrez.api_key = config.get("api_key", None)
    Entrez.sleep_between_tries = 0.34  # Respecte la limite de 3 req/s
    Entrez.max_tries = 3
    Entrez.tool = "SpookyProcessor"  # Identifie votre application

# Appeler au chargement du module
init_ncbi()

# Setup Logging
def setup_logging(output_dir=None):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    class ConsoleFilter(logging.Filter):
        def filter(self, record):
            msg = record.getMessage()
            return "processed out of" in msg or "duplicates suppressed" in msg

#    console_handler = logging.StreamHandler()
#    console_formatter = logging.Formatter("%(message)s")
#    console_handler.setFormatter(console_formatter)
#    console_handler.addFilter(ConsoleFilter())
#    console_handler.setLevel(logging.INFO)
#    logger.addHandler(console_handler)

    if output_dir:
        log_path = os.path.join(output_dir, "spooky_processor.log")
        os.makedirs(output_dir, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Masquer les logs verbeux de Biopython
    logging.getLogger("Bio").setLevel(logging.WARNING)

# Change PROTEIN_REGEX patterns to:
PROTEIN_REGEX1 = re.compile(r'^[A-Z]{2,3}_\d{5,8}(\.\d+)?$')  # Requires exactly 2-3 letters + underscore
PROTEIN_REGEX2 = re.compile(r'^[A-Z]{1,2}_\d{6,8}(\.\d+)?$')  # Requires exactly 1-2 letters + underscore
PROTEIN_REGEX3 = re.compile(r'^[A-Z]{3,6}\d{5,8}(\.\d+)?$')  # 3-6 letters followed by 5-8 digits
#ACCESSION_REGEX = re.compile(r'^([A-Z]{2,6}_\d{5,9}|[A-Z]{2}_\d{5,9})(\.\d+)?$')  # Requires exactly 2-6 letters + underscore
#ACCESSION_REGEX = re.compile(r'^([A-Z]{2,6}_\d{5,9})(\.\d+)?$')
ACCESSION_REGEX = re.compile(r'^((GCF|GCA|NZ)_[A-Z0-9]+)(\.\d+)?$')




# Constants with values from molbiotools.com
DNA_WEIGHTS = {
    'A': 313.21,
    'T': 304.19,
    'C': 289.18,
    'G': 329.21,
    'R': (313.21 + 329.21) / 2,  # A/G
    'Y': (289.18 + 304.19) / 2,  # C/T
    'S': (329.21 + 289.18) / 2,  # G/C
    'W': (313.21 + 304.19) / 2,  # A/T
    'K': (329.21 + 304.19) / 2,  # G/T
    'M': (313.21 + 289.18) / 2,  # A/C
    'B': (289.18 + 329.21 + 304.19) / 3,  # C/G/T
    'D': (313.21 + 329.21 + 304.19) / 3,  # A/G/T
    'H': (313.21 + 289.18 + 304.19) / 3,  # A/C/T
    'V': (313.21 + 289.18 + 329.21) / 3,  # A/C/G
    'N': (313.21 + 304.19 + 289.18 + 329.21) / 4  # Moyenne de A/T/C/G
}

RNA_WEIGHTS = {
    'A': 329.21,
    'U': 306.17,
    'C': 305.18,
    'G': 345.21,
    'N': (329.21 + 306.17 + 305.18 + 345.21) / 4
}

GC_CONTENT = {
    'A': 0, 'T': 0, 'C': 1, 'G': 1,
    'R': 0.5, 'Y': 0.5, 'S': 1, 'W': 0,
    'K': 0.5, 'M': 0.5, 'B': 0.67, 'D': 0.33,
    'H': 0.33, 'V': 0.67, 'N': 0.5
}

CODON_AVG_MASS = {
    # Acide aminé : Masse moyenne des codons correspondants (3nt)
    'A': 986.73,   # GCU/GCC/GCA/GCG (moyenne 4 codons)
    'C': 903.57,   # UGU/UGC
    'D': 976.70,   # GAU/GAC
    'E': 976.70,   # GAA/GAG  
    'F': 953.55,   # UUU/UUC
    'G': 986.73,   # GGU/GGC/GGA/GGG
    'H': 932.58,   # CAU/CAC
    'I': 935.56,   # AUU/AUC/AUA
    'K': 975.69,   # AAA/AAG
    'L': 935.56,   # UUA/UUG/CUU/CUC/CUA/CUG
    'M': 944.59,   # AUG
    'N': 932.58,   # AAU/AAC
    'P': 939.63,   # CCU/CCC/CCA/CCG
    'Q': 975.69,   # CAA/CAG
    'R': 1011.81,  # CGU/CGC/CGA/CGG/AGA/AGG
    'S': 903.57,   # UCU/UCC/UCA/UCG/AGU/AGC
    'T': 912.60,   # ACU/ACC/ACA/ACG
    'V': 944.59,   # GUU/GUC/GUA/GUG
    'W': 983.61,   # UGG
    'Y': 953.55,   # UAU/UAC
    '*': 0         # Stop codons
}

CODON_REL_STDDEV = {
    'A': 1.56, 'R': 0.99, 'N': 0.80, 'D': 0.80,
    'C': 0.81, 'Q': 0.84, 'E': 0.83, 'G': 1.47,
    'H': 0.83, 'I': 0.88, 'L': 1.13, 'K': 0.84,
    'M': 0.00, 'F': 0.80, 'P': 1.50, 'S': 1.34,
    'T': 1.34, 'W': 0.00, 'Y': 0.80, 'V': 1.18,
    '*': 0.80
}


ZB_GAMMA = 1e-5
FREQ_GENERATION_ERROR = 0.005  # ±0.005%
HELIX_PITCH = 3.4
HELIX_PITCH_UNCERT = 0.1
RESONANCE_COEFFICIENTS = {
    'DNA': 1.55518226281848E+17,
    'RNA': 1.61853600781306E+17,
    'mRNA': 1.50677502217234E+17,
}
MEASUREMENT_VARIANCE = (0.005) ** 2

CONFIG_FILE = get_config_path()

# NCBI Credentials
def get_ncbi_credentials(force=False):
    """Récupère les identifiants NCBI sans demander si le fichier existe"""
    config = load_ncbi_config()
    
    # 1. Vérifier d'abord si le fichier existe avec des identifiants valides
    if not force and config.get("email"):
        logging.info("Utilisation des identifiants NCBI sauvegardés")
        return config["email"], config.get("api_key")
    
    # 2. Demander seulement si nécessaire
    logging.info("Configuration NCBI requise")
    root = tk.Tk()
    root.withdraw()  # Cache la fenêtre principale vide
    
    email = simpledialog.askstring("NCBI Email", "Entrez votre email NCBI:")
    if not email:
        raise ValueError("Email requis pour utiliser les outils E-utilities")
    
    api_key = simpledialog.askstring("Clé API", "Clé API NCBI (optionnelle):", show='*')
    save_ncbi_config(email, api_key)
    return email, api_key

# Utility Functions
def clean_accession(accession):
    """Nettoie les numéros d'accession en retirant la version si présente"""
    if not accession:
        return None
    
    # Retire les espaces et convertit en majuscules
    cleaned = accession.strip().upper()
    
    # Cas des protéines sans underscore (ex: AGE13739.1 → AGE13739)
    if PROTEIN_REGEX3.match(cleaned) and '.' in cleaned:
        return cleaned.split('.')[0]
    
    """Nettoie les numéros d'accession en retirant la version si présente"""
    if not accession:
        return None
    
    # Retire les espaces et convertit en majuscules
    cleaned = accession.strip().upper()
    
    # Cas des protéines (ex: WP_123456.1 → WP_123456)
    if cleaned.count('.') == 1 and cleaned.split('.')[1].isdigit():
        return cleaned.split('.')[0]
    
    # Cas des RefSeq (ex: GCF_000001405.40 → GCF_000001405)
    if ('_' in cleaned) and cleaned.split('_')[-1].count('.') == 1:
        prefix, suffix = cleaned.split('_', 1)
        return f"{prefix}_{suffix.split('.')[0]}"
    
    return cleaned

def is_assembly_accession(acc):
    """Should be updated to explicitly include NZ_"""
    if not isinstance(acc, str):
        return False
    acc = acc.strip().upper()
#    return (acc.startswith(("GCF_", "GCA_", "NZ_")) and '_' in acc and len(acc.split('_')) == 2
    return acc.startswith(("GCF_", "GCA_", "NZ_")) and '_' in acc

def fetch_assembly_summary(accession):
    try:
        handle = Entrez.esearch(db="assembly", term=f"{accession}[Assembly Accession]", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            return None
            
        time.sleep(0.34)  # Respecter la limite de 3 requêtes/seconde
        handle = Entrez.esummary(db="assembly", id=record["IdList"][0], retmode="xml")
        summary = Entrez.read(handle)
        handle.close()
        return summary.get("DocumentSummarySet", {}).get("DocumentSummary", [{}])[0]
        
    except HTTPError as e:
        if e.code == 429:
            logging.warning("Hit NCBI rate limit, sleeping 5 seconds...")
            time.sleep(5)
            return fetch_assembly_summary(accession)  # Retry
        logging.error(f"[{accession}] Assembly summary fetch failed: {str(e)}")
        return None
    except Exception as e:
        logging.error(f"[{accession}] Assembly summary fetch failed: {str(e)}")
        return None

def download_and_extract_fasta(ftp_path, accession):
    try:
        clean_ftp_path = ftp_path.rstrip('/')
        base_name = os.path.basename(clean_ftp_path)
        fasta_url = f"https://{clean_ftp_path.replace('ftp://', '')}/{base_name}_genomic.fna.gz"
        
        response = requests.get(fasta_url, stream=True, timeout=60)
        response.raise_for_status()
        
        with gzip.GzipFile(fileobj=BytesIO(response.content)) as gz_file:
            content = gz_file.read().decode('utf-8')
        
        records = list(SeqIO.parse(StringIO(content), "fasta"))
        sequence = "".join(str(r.seq) for r in records).upper()
        description = records[0].description
        
        logging.info(f"[{accession}] Successfully downloaded assembly")
        # Retourner 4 valeurs comme attendu par fetch_sequence()
        return sequence, description, "DNA", False  # Ajout du type et circular=False
    
    except Exception as e:
        logging.error(f"[{accession}] FASTA download failed: {str(e)}")
        logging.debug(f"Used URL: {fasta_url}", exc_info=True)
        return None, None, "DNA", False  # Même en cas d'erreur, retourner 4 valeurs

# Scientific Calculations
def calculate_gc_percent(sequence):
    """Calcul du %GC en tenant compte des bases ambiguës"""
    # Remplacer 'GC_CONTENT' par 'DNA_WEIGHTS' pour la validation des bases
    total_gc = sum(GC_CONTENT.get(base, 0.5) for base in sequence.upper() if base in DNA_WEIGHTS)
    total = sum(1 for base in sequence.upper() if base in DNA_WEIGHTS)
    return round((total_gc / total * 100), 2) if total > 0 else 0

def count_valid_bases(sequence):
    """Compte uniquement les bases valides"""
    return sum(1 for base in sequence.upper() if base in DNA_WEIGHTS)

def calculate_molecular_weight(sequence, is_rna=False):
    weight_dict = RNA_WEIGHTS if is_rna else DNA_WEIGHTS
    seq_upper = sequence.upper()
    total_mass = sum(weight_dict.get(base, weight_dict["N"]) for base in seq_upper)
    
    # Version avec math au lieu de numpy
    mass_abs_uncert = math.sqrt(total_mass * MEASUREMENT_VARIANCE)
    mass_rel_uncert = (mass_abs_uncert / total_mass) * 100 if total_mass != 0 else 0
    zb_uncert_percent = ZB_GAMMA * math.log10(total_mass) * 100
    mass_total_rel_uncert = math.sqrt(mass_rel_uncert**2 + FREQ_GENERATION_ERROR**2 + zb_uncert_percent**2)
    da_per_nt = total_mass / len(sequence) if sequence else 0
    return round(total_mass, 2), round(da_per_nt, 5), round(mass_total_rel_uncert, 2)

def calculate_resonance_frequency(nt, seq_type="DNA"):
    coefficient = RESONANCE_COEFFICIENTS.get(seq_type, RESONANCE_COEFFICIENTS['DNA'])
    freq_res = coefficient / nt
    rel_abs_uncertainty = np.sqrt((HELIX_PITCH_UNCERT / HELIX_PITCH) ** 2)
    res_total_rel_uncert = np.sqrt(rel_abs_uncertainty ** 2 + FREQ_GENERATION_ERROR ** 2) * 100
    return round(freq_res, 2), round(res_total_rel_uncert, 2)

# Sequence fetching
def resolve_refseq_link(accession):
    try:
        handle = Entrez.elink(dbfrom="nucleotide", id=accession, linkname="nuccore")
        results = Entrez.read(handle)
        handle.close()
        if results and len(results[0]["LinkSetDb"]) > 0:
            linked_id = results[0]["LinkSetDb"][0]["Link"][0]["Id"]
            return linked_id
    except Exception as e:
        pass
    try:
        handle = Entrez.esummary(db="nucleotide", id=accession)
        summary = Entrez.read(handle)[0]
        handle.close()
        if "Comment" in summary:
            match = re.search(r"\b(NM_\d+\.\d+|NC_\d+\.\d+|NZ_.+)\b", summary["Comment"])
            if match:
                return match.group(1)
    except Exception as e:
        pass
    return accession

local_cache = {}

def cached_efetch(db, id, **kwargs):
    cache_key = f"{db}_{id}"
    if cache_key in local_cache:
        return local_cache[cache_key]
    
    try:
        handle = Entrez.efetch(db=db, id=id, **kwargs)
        result = handle.read()
        handle.close()
        local_cache[cache_key] = result
        return result
    except Exception as e:
        logging.error(f"Cached efetch failed: {str(e)}")
        raise

def get_sequence_type(accession):
    """Determine sequence type with strict validation"""
    if not isinstance(accession, str):
        return None

    accession = accession.strip().upper()

    if not accession or '_' not in accession:
        # Check for protein accessions without underscore
        if PROTEIN_REGEX3.match(accession):
            return "protein"
        return None

    if not isinstance(accession, str):
        return None

    accession = accession.strip().upper()

    if not accession or '_' not in accession:
        return None

    # Refuse explicitement les faux formats GC12345 et GCA12345
    if accession.startswith("GC") and '_' not in accession:
        return None

    # 1. ARN d'abord
    if accession.startswith(('NR_', 'XR_', 'NM_', 'XM_')):
        return "rna"

    # 2. Assemblages
    if is_assembly_accession(accession):
        return "dna"

    # 3. ADN : traiter explicitement GC_
    if accession.startswith('GC_'):
        return "dna"

    # 4. ADN générique
    if ACCESSION_REGEX.match(accession):
        return "dna"

    # 5. Protéines (en dernier)
    protein_prefixes = ('NP_', 'XP_', 'YP_', 'WP_', 'AP_')
    if accession.startswith(protein_prefixes) or PROTEIN_REGEX1.match(accession) or PROTEIN_REGEX2.match(accession):
        return "protein"

    return None

def fetch_sequence(accession):
    """Version optimisée avec nettoyage systématique de l'accession"""
    # ==============================================
    # 0. NETTOYAGE DE L'ACCESSION
    # ==============================================
    original_accession = accession
    accession = clean_accession(accession)
    if not accession:
        logging.error(f"[{original_accession}] Invalid accession after cleaning")
        return None, "", "DNA", False
    
    # ==============================================
    # 1. VALIDATION ET DÉTECTION DU TYPE
    # ==============================================
    seq_type = get_sequence_type(accession)
    if not seq_type:
        logging.error(f"[{accession}] Invalid accession format")
        return None, "", "DNA", False

    # ==============================================
    # 2. TRAITEMENT PROTÉINE (ancien is_protein)
    # ==============================================
    if seq_type == "protein":
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # A. Récupération séquence AA
                handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
                protein_record = SeqIO.read(handle, "fasta")
                handle.close()

                # B. Tentative de trouver le CDS parent
                try:
                    clean_id = accession.split('.')[0]
                    handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=clean_id, linkname="protein_nuccore")
                    links = Entrez.read(handle)
                    handle.close()
                    
                    if links[0]["LinkSetDb"]:
                        nucl_id = links[0]["LinkSetDb"][0]["Link"][0]["Id"]
                        handle = Entrez.efetch(db="nucleotide", id=nucl_id, rettype="gb", retmode="text")
                        nucl_record = SeqIO.read(handle, "genbank")
                        handle.close()
                        
                        for feature in nucl_record.features:
                            if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                                if clean_id in feature.qualifiers["protein_id"][0]:
                                    cds_seq = str(feature.extract(nucl_record.seq))
                                    return cds_seq, f"CDS for {accession}", "DNA", False
                except Exception as cds_error:
                    logging.warning(f"[{accession}] CDS fetch failed, using protein seq: {cds_error}")

                # C. Fallback sur la séquence AA
                return str(protein_record.seq), protein_record.description, "protein", False

            except HTTPError as e:
                if e.code == 400 and attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                logging.error(f"[{accession}] Protein HTTP error {e.code}")
                return None, "", "protein", False
            except Exception as e:
                logging.error(f"[{accession}] Protein fetch error: {str(e)}")
                return None, "", "protein", False

    # ==============================================
    # 3. CAS DES GENE IDs (chiffres uniquement) - SECTION COMPLÉTÉE
    # ==============================================
    elif accession.isdigit():
        try:
            # Récupérer les informations du gène
            handle = Entrez.efetch(db="gene", id=accession, rettype="xml")
            gene_info = Entrez.read(handle)[0]
            handle.close()
            
            # Essayer de trouver le nucléotide lié
            if "Entrezgene_locus" in gene_info and gene_info["Entrezgene_locus"]:
                for locus in gene_info["Entrezgene_locus"]:
                    if "Gene-commentary_accession" in locus:
                        nucl_accession = locus["Gene-commentary_accession"]
                        logging.info(f"[{accession}] Found linked nucleotide: {nucl_accession}")
                        
                        # Récupérer la séquence nucléotide
                        handle = Entrez.efetch(db="nucleotide", id=nucl_accession, 
                                             rettype="fasta", retmode="text")
                        record = SeqIO.read(handle, "fasta")
                        handle.close()
                        
                        return str(record.seq).upper(), record.description, "DNA", False
            
            # Fallback: essayer de trouver via elink
            handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=accession)
            links = Entrez.read(handle)
            handle.close()
            
            if links and links[0]["LinkSetDb"]:
                nucl_id = links[0]["LinkSetDb"][0]["Link"][0]["Id"]
                handle = Entrez.efetch(db="nucleotide", id=nucl_id, 
                                     rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(record.seq).upper(), record.description, "DNA", False
                
            logging.error(f"[{accession}] No linked nucleotide found for gene")
            return None, "", "DNA", False

        except HTTPError as e:
            logging.error(f"[{accession}] Gene HTTP error {e.code}")
            return None, "", "DNA", False
        except Exception as e:
            logging.error(f"[{accession}] Gene fetch error: {str(e)}")
            return None, "", "DNA", False

    # ==============================================
    # 4. CAS DES GÉNOMES (GCF_/GCA_/NZ_)
    # ==============================================
    elif is_assembly_accession(accession):
        try:
            summary = fetch_assembly_summary(accession)
            if summary:
                ftp_path = summary.get("FtpPath_RefSeq", "") or summary.get("FtpPath_GenBank", "")
                if ftp_path:
                    sequence, description, seq_type, is_circular = download_and_extract_fasta(ftp_path, accession)
                    return sequence, description, seq_type, is_circular
        except Exception as e:
            logging.error(f"[{accession}] Assembly fetch failed: {str(e)}")
        return None, "", "DNA", False

    # ==============================================
    # 5. CAS PAR DÉFAUT (nucléotides standards)
    # ==============================================
    else:
        max_retries = 3
        for attempt in range(max_retries):
            try:
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(record.seq).upper(), record.description, "DNA", "circular" in record.description.lower()
            except HTTPError as e:
                if attempt < max_retries - 1:
                    time.sleep((2 ** attempt) * 0.5)
                    continue
                logging.error(f"[{accession}] Nucleotide HTTP error {e.code}")
                return None, "", "DNA", False
            except Exception as e:
                logging.error(f"[{accession}] Nucleotide fetch error: {str(e)}")
                return None, "", "DNA", False


    logging.error(f"[{accession}] All fetch methods failed")
    return None, "", "DNA", False

def fetch_protein_sequence(protein_id):
    """Récupère la séquence protéique depuis NCBI."""
    try:
        handle = Entrez.efetch(db="protein", id=protein_id.split('.')[0], rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        logging.error(f"[{protein_id}] Protein fetch failed: {e}")
        return None

def calculate_from_cds(cds_seq, description):
    """Calcule la masse moléculaire depuis une CDS."""
    total_mass, da_per_nt, _ = calculate_molecular_weight(cds_seq)
    return {
        "description": description[:200],
        "mass": int(total_mass),
        "nt": len(cds_seq),
        "da_per_nt": round(da_per_nt, 1),
        "is_circular": False
    }

def calculate_from_protein(aa_seq, description):
    """Fallback : calcule la masse depuis la séquence AA et ajoute l’incertitude relative."""
    avg_mass = sum(CODON_AVG_MASS.get(aa, 0) for aa in aa_seq)
    estimated_nt = len(aa_seq) * 3

    # Calculer l’incertitude relative combinée
    rel_variances = [
        (CODON_REL_STDDEV.get(aa, 0) / 100)**2
        for aa in aa_seq
    ]
    total_rel_variance = sum(rel_variances)
    rel_uncertainty = (total_rel_variance ** 0.5) * 100  # en %

    return {
        "description": description[:200],
        "mass": int(avg_mass),
        "nt": estimated_nt,
        "da_per_nt": round(avg_mass / estimated_nt, 1),
        "is_circular": False,
        "relative_uncertainty": round(rel_uncertainty, 3)  # % final
    }

def fetch_protein_data(protein_id):
    """Orchestre le tout: essaie CDS, sinon fallback."""
    # 1. Récupérer la protéine
    protein_record = fetch_protein_sequence(protein_id)
    if not protein_record:
        return None

    # 2. Essayer CDS (sauf pour WP_/XP_ si optimisation souhaitée)
    if not protein_id.startswith(('WP_', 'XP_')):
        cds_seq, cds_desc = fetch_cds_for_protein(protein_id)
        if cds_seq:
            return calculate_from_cds(cds_seq, cds_desc)

    # 3. Fallback sur la protéine
    return calculate_from_protein(str(protein_record.seq), protein_record.description)

            
def fetch_cds_for_protein(protein_id):
    """Tente de trouver la CDS correspondant à une protéine."""
    try:
        clean_id = protein_id.split('.')[0]
        handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=clean_id, linkname="protein_nuccore")
        links = Entrez.read(handle)
        handle.close()
        if links[0]["LinkSetDb"]:
            nucl_id = links[0]["LinkSetDb"][0]["Link"][0]["Id"]
            handle = Entrez.efetch(db="nucleotide", id=nucl_id, rettype="gb", retmode="text")
            nucl_record = SeqIO.read(handle, "genbank")
            handle.close()
            # Cherche la CDS correspondante
            for feature in nucl_record.features:
                if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                    if clean_id in feature.qualifiers["protein_id"][0]:
                        cds_seq = str(feature.extract(nucl_record.seq))
                        if not cds_seq:
                            logging.info(f"[{protein_id}] No CDS sequence found, using protein fallback")
                            return None, ""
                        return cds_seq, f"CDS for {protein_id}"
    except Exception as e:
        logging.warning(f"[{protein_id}] CDS fetch failed: {e}")
    return None, ""


def process_accession(row, args):
    """Version optimisée avec :
    - Gestion robuste des accessions
    - Validation centralisée
    - Gestion modulaire des types
    - Logging amélioré
    """

    raw_accession = row["Assembly Accession"].strip()
    accession = clean_accession(raw_accession)
    
    if not accession:
        logging.error(f"Invalid accession: {raw_accession}")
        return None, None
        
    logging.debug(f"Processing {accession} (original: {raw_accession})")  
    accession = row["Assembly Accession"].strip()
    
    # Validation
    if not accession:
        logging.error("Empty accession")
        return None, None
        
    seq_type = get_sequence_type(accession)
    logging.debug(f"[{accession}] Detected type: {seq_type}")

    # Protéines
    if seq_type == "protein":
        result = fetch_protein_data(accession)
        if not result:
            return None, None
            
        return format_protein_result(result, accession, args)

    # ADN/ARN
    return process_dna_rna(accession, args)

def format_protein_result(result, accession, args):
    line_mw = [
        f'PROTEIN {accession} | {result["description"]}',
        "CUST",                                             
        "",
        f'Da/nt={result["mass"]/result["nt"]:.5f} | nt={result["nt"]} | Uncertainty=±{result.get("relative_uncertainty")}%',
        f"M{result['mass']}",                                # ✅ sans guillemets
        "", "", 180                                           # ✅ sans guillemets
    ]

    # Ligne RES (si activée)
    line_res = None
    if args.get("include_resonance"):
        line_res = [
            f'PROTEIN {accession}',
            "CUST",                                            # ✅ sans guillemets
            "",
            f'Uncertainty ±{result.get("relative_uncertainty", 0):.2f}%',
            f"BP{result['nt']}",                               # ✅ sans guillemets
            "", "", 180                                         # ✅ sans guillemets
        ]
    
    return line_mw, line_res


def process_dna_rna(accession, args):
    sequence, desc, seq_type, is_circular = fetch_sequence(accession)
    if not sequence:
        return None, None
    
    # Calcul des propriétés
    nt = len(sequence)
    gc_percent = calculate_gc_percent(sequence)
    total_mass, da_per_nt, mass_uncert = calculate_molecular_weight(sequence, seq_type == "rna")
    n_count = sequence.upper().count('N')
    n_percent = (n_count / nt) * 100 if nt > 0 else 0
    _, res_uncert = calculate_resonance_frequency(nt, seq_type)

    # Ligne MW
    line_mw = [
        f'{seq_type.upper()} {accession} | {desc}',
        "CUST",     
        "",
        f'GC%={gc_percent:.2f}% | N={n_count} ({n_percent:.2f}%) | Da/nt={da_per_nt:.5f} | nt={nt} | MW uncertainty ±{mass_uncert:.2f}%',
        f"M{total_mass:.1f}", 
        "", "", 180               
    ]

    # Ligne RES
    line_res = None
    if args.get("include_resonance"):
        prefix_map = {
            ('DNA', True): 'BC',
            ('DNA', False): 'BL',
            ('RNA', True): 'BCR', 
            ('RNA', False): 'BLR',
            ('mRNA', True): 'BCm',
            ('mRNA', False): 'BLm'
        }
        prefix = prefix_map.get((seq_type, is_circular), 'BXX')
    
        line_res = [
            f'{desc} ({seq_type} {"circular" if is_circular else "linear"})',  # ✅ guillemets
            "CUST",                                                              # ✅ sans guillemets
            "",
            f'Uncertainty ±{res_uncert:.2f}%',                                 # ✅ guillemets
            f"{prefix}{len(sequence)}",                                          # ✅ sans guillemets
            "", "", 180                                                          # ✅ sans guillemets
        ]  
        return line_mw, line_res
   
def read_tsv_file(file_path, include_gca=False):
    """Lit un fichier TSV ou CSV et retourne une liste d'accessions valides."""
    accessions = []
    try:
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            if not lines:
                logging.error("Fichier vide")
                return []

            first_line = lines[0].lower()
            has_header = any(keyword in first_line for keyword in ['assembly', 'accession', 'genome'])
            delimiter = '\t' if '\t' in first_line else ','

            acc_col_index = 0
            if has_header:
                headers = [h.strip().lower() for h in first_line.split(delimiter)]
                for i, header in enumerate(headers):
                    if 'assembly' in header or 'accession' in header:
                        acc_col_index = i
                        break
                data_lines = lines[1:]
            else:
                data_lines = lines

            unique_accessions = OrderedDict()
            for line in data_lines:
                if line.startswith('#'):
                    continue
                parts = [p.strip() for p in line.split(delimiter)]
                if len(parts) > acc_col_index:
                    acc = parts[acc_col_index]
                    if acc and acc.lower() != "nan":
                        cleaned = clean_accession(acc)
                        if cleaned and (include_gca or not cleaned.startswith("GCA_")):
                            unique_accessions[cleaned] = True
            return [{"Assembly Accession": acc} for acc in unique_accessions.keys()]
    except Exception as e:
        logging.critical(f"Erreur lors de la lecture du fichier: {str(e)}")
        return []

def save_output(file_path, lines):
    """Sauvegarde les résultats sans échapper les guillemets"""
    try:
        with open(file_path, "w", newline='', encoding="utf-8") as f:
            writer = csv.writer(f, 
                              quoting=csv.QUOTE_MINIMAL,  # Quote only when needed
                              escapechar='',             # No escape character
                              doublequote=False)         # Don't double quotes
            
            for line in lines:
                # Formatage manuel des champs nécessitant des guillemets
                formatted = [
                    line[0] if ',' not in line[0] else f'"{line[0]}"',  # Description
                    line[1],  # CUST
                    '',       # Champ vide
                    line[3] if ',' not in line[3] else f'"{line[3]}"',  # Metadata
                    line[4],  # M/Bxxx
                    '', '',  # Champs vides
                    line[7]  # 180
                ]
                writer.writerow(formatted)
    except Exception as e:
        logging.critical(f"Failed to save {file_path}: {str(e)}")
        raise

from concurrent.futures import TimeoutError

def process_accession_wrapper(row, args):
    """Wrapper avec timeout et gestion d'erreur fine"""
    try:
        # Timeout pour chaque requête NCBI
        Entrez.timeout = 30  # 30s max par requête NCBI
        
        # Exécution
        line_mw, line_res = process_accession(row, args)
        return line_mw, line_res, row["Assembly Accession"]
    except Exception as e:
        logging.error(f"[{row['Assembly Accession']}] Wrapper error: {str(e)}")
        return None, None, row["Assembly Accession"]

def timeout_handler(signum, frame):
    """Handler pour le timeout"""
    raise TimeoutError("Processing timeout after 5 minutes")

def process_tsv(input_file, output_dir, include_gca=False, max_workers=4, include_resonance=False):
    """Version finale avec gestion robuste des threads et timeout global"""
    # Configuration initiale
    failed_accessions = []
    executor = None
    
    try:
        # 1. Initialisation avec timeout
        socket.setdefaulttimeout(60)  # Timeout réseau global (60s)
        setup_logging(output_dir)
        
        # 2. Préparation des fichiers de sortie
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        mw_output_path = os.path.join(output_dir, f"{base_name}_spooky_mw_{timestamp}.csv")
        res_output_path = os.path.join(output_dir, f"{base_name}_spooky_res_{timestamp}.csv") if include_resonance else None
        
        os.makedirs(output_dir, exist_ok=True)

        # 3. Lecture du fichier d'entrée
        df = read_tsv_file(input_file, include_gca)
        if not df:
            logging.error("No valid accessions found")
            return False

        # 4. Configuration du ThreadPool avec gestionnaire de contexte
        with ThreadPoolExecutor(
            max_workers=min(max_workers, 8),  # Limite à 8 threads max
            thread_name_prefix="spooky_worker"
        ) as executor:
            
            # 5. Soumission des tâches avec timeout individuel
            future_to_acc = {
                executor.submit(
                    process_accession_wrapper, 
                    row, 
                    {"include_resonance": include_resonance}
                ): row["Assembly Accession"] 
                for row in df
            }

            results = []
            seen_accessions = set()
            
            for future in as_completed(future_to_acc):
                acc = future_to_acc[future]
                try:
                    line_mw, line_res, _ = future.result(timeout=300)  # 5 min par accession
                    if line_mw and acc not in seen_accessions:
                        seen_accessions.add(acc)
                        results.append((line_mw, line_res))
                except TimeoutError:
                    logging.error(f"[{acc}] Timeout after 5 minutes")
                    failed_accessions.append(acc)
                except Exception as e:
                    logging.error(f"[{acc}] Error: {str(e)}")
                    failed_accessions.append(acc)

        # 6. Dédoublonnage et sauvegarde (hors du with pour libérer les threads)
        if results:
            save_results(mw_output_path, res_output_path, results, include_resonance)

        # 7. Rapport final
        log_final_report(len(results), len(df), failed_accessions)
        return True

    except Exception as e:
        logging.critical(f"Fatal error: {str(e)}", exc_info=True)
        return False
    finally:
        # Nettoyage garantie
        if executor:
            executor.shutdown(wait=False)  # Arrêt immédiat des threads
        import gc
        gc.collect()


# Fonctions auxiliaires séparées pour plus de clarté
def save_results(mw_path, res_path, results, include_resonance):
    """Sauvegarde les résultats dans les fichiers"""
    unique_masses = {}
    unique_res = {}
    
    for line_mw, line_res in results:
        # Dédoublonnage masses
        mass_key = round(float(line_mw[4][1:]), 1)  # Extrait M{value}
        if mass_key not in unique_masses:
            unique_masses[mass_key] = line_mw
        
        # Dédoublonnage résonances
        if include_resonance and line_res:
            res_key = line_res[4]  # BL/BC code
            if res_key not in unique_res:
                unique_res[res_key] = line_res
    
    # Sauvegarde
    with open(mw_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        writer.writerows(unique_masses.values())
    
    if include_resonance and unique_res:
        with open(res_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerows(unique_res.values())

def log_final_report(success_count, total_count, failed_accessions):
    """Journalise le rapport final"""
    logging.info(f"Successfully processed {success_count}/{total_count} accessions")
    if failed_accessions:
        logging.warning(f"Failed accessions ({len(failed_accessions)}): {', '.join(failed_accessions)}")

# CLI Entry Point
def parse_args():
    parser = argparse.ArgumentParser(description="Process NCBI accessions")
    parser.add_argument("tsv_file", help="Input TSV file path")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--include-gca", action="store_true", help="Include GCA_ accessions")
    parser.add_argument("--include-resonance", action="store_true", help="Generate resonance frequency CSV")
    parser.add_argument("--max-workers", type=int, default=4, help="Max parallel threads")
    parser.add_argument("--reset-ncbi", action="store_true", help="Force re-entry of NCBI credentials")
    return parser.parse_args()

# GUI
class SpookyProcessorGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Spooky Processor v2.0")
        self.geometry("900x700")
        self.resizable(True, True)
        self.setup_ui()
        self.setup_logging()
        self.show_readme()

    def setup_ui(self):
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Input Section
        input_frame = ttk.LabelFrame(main_frame, text="Input", padding="10")
        input_frame.pack(fill=tk.X, pady=5)
        ttk.Label(input_frame, text="TSV File:").grid(row=0, column=0, sticky=tk.W)
        self.input_entry = ttk.Entry(input_frame, width=60)
        self.input_entry.grid(row=0, column=1, padx=5)
        ttk.Button(input_frame, text="Browse", command=self.select_tsv).grid(row=0, column=2)

        # Output Section
        output_frame = ttk.LabelFrame(main_frame, text="Output", padding="10")
        output_frame.pack(fill=tk.X, pady=5)
        ttk.Label(output_frame, text="Directory:").grid(row=0, column=0, sticky=tk.W)
        self.output_entry = ttk.Entry(output_frame, width=60)
        self.output_entry.grid(row=0, column=1, padx=5)
        ttk.Button(output_frame, text="Browse", command=self.select_out).grid(row=0, column=2)

        # Options
        options_frame = ttk.LabelFrame(main_frame, text="Options", padding="10")
        options_frame.pack(fill=tk.X, pady=5)
        self.include_gca = tk.BooleanVar()
        self.include_resonance = tk.BooleanVar(value=True)
        self.max_workers = tk.IntVar(value=4)
        ttk.Checkbutton(options_frame, text="Include GCA_", variable=self.include_gca).grid(row=0, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Include Resonance Data", variable=self.include_resonance).grid(row=0, column=1, sticky=tk.W)
        ttk.Label(options_frame, text="Threads:").grid(row=1, column=0, sticky=tk.W)
        ttk.Spinbox(options_frame, from_=1, to=16, textvariable=self.max_workers, width=4).grid(row=1, column=1, sticky=tk.W)

        # Log Area
        log_frame = ttk.LabelFrame(main_frame, text="Processing Log", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.log_text = tk.Text(log_frame, height=15, state=tk.DISABLED)
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Progress & Controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=5)
        self.progress = ttk.Progressbar(control_frame, mode='indeterminate')
        self.progress.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        ttk.Button(control_frame, text="Start", command=self.start_processing).pack(side=tk.RIGHT)

    def setup_logging(self):
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()

        # Handler fichier
        log_path = os.path.join(os.path.dirname(__file__), "spooky_processor.log")
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)

        # Handler GUI
        class GuiLogHandler(logging.Handler):
            def __init__(self, text_widget):
                super().__init__()
                self.text_widget = text_widget

            def emit(self, record):
                with log_lock:
                    msg = self.format(record)
                    self.text_widget.config(state=tk.NORMAL)
                    self.text_widget.insert(tk.END, msg + "\n")
                    self.text_widget.see(tk.END)
                    self.text_widget.config(state=tk.DISABLED)
        self.gui_handler = GuiLogHandler(self.log_text)
        self.gui_handler.setLevel(logging.INFO)
        self.gui_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(self.gui_handler)

    def select_tsv(self):
        path = filedialog.askopenfilename(filetypes=[("TSV files", "*.tsv")])
        if path:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, path)

    def select_out(self):
        path = filedialog.askdirectory()
        if path:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, path)

    def show_readme(self):
        readme_path = os.path.join(os.path.dirname(__file__), "README.txt")
        if os.path.exists(readme_path):
            try:
                with open(readme_path, "r", encoding="utf-8") as f:
                    content = f.read()
                self.log_text.config(state=tk.NORMAL)
                self.log_text.insert(tk.END, f"README:\n{content}\n\n")
                self.log_text.config(state=tk.DISABLED)
            except Exception as e:
                logging.warning(f"Could not read README: {str(e)}")

    def get_credentials(self):
        config = load_ncbi_config()
        if config.get("email") and not self.include_gca.get():
            if messagebox.askyesno("Credentials", "Use saved credentials?"):
                return config["email"], config.get("api_key")
        email = simpledialog.askstring("NCBI Email", "Enter your NCBI email:")
        if not email:
            return None, None
        api_key = simpledialog.askstring("API Key", "Optional API key:", show='*')
        save_ncbi_config(email, api_key)
        return email, api_key

    def start_processing(self):
        try:
            input_file = self.input_entry.get()
            output_dir = self.output_entry.get()
        
            if not all([input_file, output_dir]):
                messagebox.showerror("Error", "Please specify input and output paths.")
                return

            # Configuration NCBI avant traitement (surcharge les valeurs par défaut)
            try:
                email, api_key = get_ncbi_credentials()
                Entrez.email = email
                if api_key:
                    Entrez.api_key = api_key
                Entrez.sleep_between_tries = 1  # Plus agressif en GUI
            except Exception as e:
                messagebox.showerror("Error", f"NCBI config failed: {str(e)}")
                return

            self.progress.start()
            threading.Thread(
                target=self.run_processing,
                args=(input_file, output_dir),
                daemon=True
            ).start()
        
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {str(e)}")
            self.progress.stop()

    def run_processing(self, input_file, output_dir):
        try:
            # Ajouter un timeout pour les requêtes NCBI
            Entrez.timeout = 30  # 30 secondes max par requête
        
            process_tsv(
                input_file,
                output_dir,
                include_gca=self.include_gca.get(),
                max_workers=self.max_workers.get(),
                include_resonance=self.include_resonance.get()
            )
        except Exception as e:
            self.after(100, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.progress.stop()
            self.after(100, lambda: messagebox.showinfo("Info", "Traitement terminé"))
            # Force le nettoyage
            import gc
            gc.collect()
# Point d'entrée principal
# Deux modes d'exécution :
# Mode normal : python spooky_processor.py (lance la GUI)
# Mode test : python spooky_processor.py --test (lance les tests)

def main():
    """Point d'entrée principal pour l'interface graphique"""
    # Configuration initiale d'Entrez
    config_ncbi()
    app = SpookyProcessorGUI()
    app.mainloop()

def config_ncbi():
    """Configure les paramètres globaux d'Entrez"""
    Entrez.sleep_between_tries = 15  # Secondes entre les tentatives
    Entrez.max_tries = 3  # Nombre maximum de tentatives
    Entrez.email = "your_email@example.com"  # Doit être une adresse valide
    
    # Charger la clé API depuis la config si disponible
    config = load_ncbi_config()
    if config.get("api_key"):
        Entrez.api_key = config["api_key"]

def run_tests():
    """Exécute tous les tests unitaires avec les corrections nécessaires"""
    test_results = {
        'passed': 0,
        'failed': 0,
        'errors': 0
    }

    def assert_test(input_val, expected, test_func, description):
        """Fonction helper pour exécuter un test"""
        nonlocal test_results
        try:
            result = test_func(input_val)
            if result == expected:
                test_results['passed'] += 1
                print(f"✅ {description}")
                return True
            else:
                test_results['failed'] += 1
                print(f"❌ {description} (obtenu: {result}, attendu: {expected})")
                return False
        except Exception as e:
            test_results['errors'] += 1
            print(f"⚠️ ERREUR sur {description}: {str(e)}")
            return False

    print("\n=== TESTS UNITAIRES ===")
    
    # Test clean_accession()
    print("\n=== TEST CLEAN_ACCESSION ===")
    clean_accession_cases = [
        ("GCF_000001405.40", "GCF_000001405", "RefSeq avec version"),
        ("WP_123456.1", "WP_123456", "Protéine standard"),
        ("NC_123456.789", "NC_123456", "Chromosome"),
        ("NZ_ABCD123456.1", "NZ_ABCD123456", "Génome assembleur"),
        ("GCA_123.456", "GCA_123", "GenBank avec version"),
        ("INVALID", "INVALID", "Chaîne invalide"),
        ("123.456", "123", "Numérique avec décimale"),
        ("", None, "Chaîne vide"),
        (None, None, "Valeur None")
    ]
    
    for input_acc, expected, desc in clean_accession_cases:
        assert_test(input_acc, expected, clean_accession, desc)

    # Test get_sequence_type()
    print("\n=== TEST GET_SEQUENCE_TYPE ===")
    sequence_type_cases = [
        ("GCF_000001405", "dna", "RefSeq standard"),
        ("GCA_000001405", "dna", "GenBank standard"),
        ("NZ_ABCD123456", "dna", "Génome assembleur"),
        ("WP_123456.1", "protein", "Protéine avec version"),
        ("NP_123456", "protein", "Protéine RefSeq"),
        ("NM_123456", "rna", "ARN messager"),
        ("NR_123456", "rna", "ARN non-codant"),
        ("GC_12345", "dna", "Format GC_ valide"),
        ("GCA_12345", "dna", "Format GCA_ valide"),
        ("GC12345", None, "Format GC sans underscore"),
        ("GCA12345", None, "Format GCA sans underscore"),
        ("INVALID_ACC", None, "Accession invalide"),
        ("RANDOMTEXT", None, "Texte aléatoire"),
        ("123456789", None, "Numérique seul"),
        ("", None, "Chaîne vide"),
        (None, None, "Valeur None"),
        (12345, None, "Entier"),
        (["GCF_123"], None, "Liste (type invalide)")
    ]
    
    for acc, expected, desc in sequence_type_cases:
        assert_test(acc, expected, get_sequence_type, desc)

    # Section 3: Test read_tsv_file()
    print("\n=== TEST READ_TSV_FILE ===")
    test_files = {
        "with_header.tsv": "Assembly Accession\nGCF_001.1\nGCF_002.2\nINVALID",
        "no_header.tsv": "GCF_003.3\nGCF_004.4\n\n# Comment",
        "empty_file.tsv": "",
        "mixed_columns.tsv": "ID\tAssembly\n1\tGCF_005.5\n2\tGCF_006.6"
    }
    
    for filename, content in test_files.items():
        try:
            with open(filename, 'w') as f:
                f.write(content)
            
            print(f"\nTesting file: {filename}")
            results = read_tsv_file(filename)
            
            if not results:
                print("  ⚠ No results found")
                continue
                
            for i, row in enumerate(results, 1):
                acc = row.get("Assembly Accession", "MISSING_KEY")
                print(f"  {i}. {acc} (Cleaned: {clean_accession(acc)})")
                
        except Exception as e:
            print(f"  ❌ Error processing {filename}: {str(e)}")

    print(f"\nRÉSUMÉ: {test_results['passed']} succès, "
          f"{test_results['failed']} échecs, "
          f"{test_results['errors']} erreurs")
    return test_results['failed'] == 0 and test_results['errors'] == 0

if __name__ == "__main__":
    import sys
    import os  # ✅ importer os ici

    if '--test' in sys.argv:
        print("=== DÉBUT DES TESTS ===")
        if not run_tests():
            os._exit(1)  # 🔥 tue immédiatement
        else:
            os._exit(0)
    else:
        main()

    print("✅ Fin du script (forcé)", flush=True)
    os._exit(0)  # ✅ tue le processus après main()


