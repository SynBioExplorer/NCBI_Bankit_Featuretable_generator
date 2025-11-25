#!/usr/bin/env python3
"""
Convert GFF3 file to NCBI BankIt Feature Table format.

NCBI Feature Table format:
- Tab-delimited, 5 columns
- Column 1: Start location
- Column 2: Stop location
- Column 3: Feature name
- Column 4: Qualifier name
- Column 5: Qualifier value
"""

import re
import sys
from collections import defaultdict


# Features to include in NCBI submission
NCBI_FEATURES = {
    'gene', 'CDS', 'mRNA', 'tRNA', 'rRNA', 'ncRNA',
    'misc_feature', 'rep_origin', 'repeat_region',
    'exon', 'intron', 'promoter', 'terminator',
    'regulatory', 'misc_binding', 'primer_bind', 'primer_bind_reverse',
    'ARS', 'ARS_consensus_sequence', 'modified_base'
}

# Map non-standard feature names to NCBI-accepted names
FEATURE_NAME_MAP = {
    'primer_bind_reverse': 'primer_bind',
    'ARS': 'rep_origin',
    'ARS_consensus_sequence': 'misc_feature',
}


def parse_attributes(attr_string):
    """Parse GFF3 attributes column into a dictionary."""
    attrs = {}
    if not attr_string or attr_string == '.':
        return attrs

    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            # URL decode common characters
            value = value.replace('%2C', ',').replace('%3B', ';').replace('%2F', '/')
            attrs[key] = value
    return attrs


def extract_gene_name(name):
    """Extract gene name from feature Name attribute.

    Examples:
        'PFS2 mRNA' -> 'PFS2'
        'PFS2 CDS' -> 'PFS2'
        'PFS2 gene' -> 'PFS2'
        'tRNA-Phe' -> 'tRNA-Phe'
    """
    if not name:
        return None

    # Remove common suffixes
    for suffix in [' mRNA', ' CDS', ' gene', ' ncRNA', ' tRNA', ' rRNA']:
        if name.endswith(suffix):
            return name[:-len(suffix)]

    return name


def parse_gff(gff_file):
    """Parse GFF3 file and return features grouped by ID for multi-interval features."""
    features = []
    multi_interval_features = defaultdict(list)  # Group by ID for spliced features
    intron_regions = set()  # Track intron coordinates to exclude from CDS

    # First pass: collect all intron regions
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            if feature_type == 'intron':
                intron_regions.add((int(start), int(end)))

    # Second pass: parse all features
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Skip features not relevant for NCBI
            if feature_type not in NCBI_FEATURES:
                continue

            attrs = parse_attributes(attributes)
            feature_id = attrs.get('ID')
            name = attrs.get('Name', '')

            feature = {
                'seqid': seqid,
                'type': feature_type,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'phase': phase,
                'name': name,
                'id': feature_id,
                'attrs': attrs
            }

            # Group multi-interval features by ID
            if feature_id and '.' in feature_id:
                multi_interval_features[feature_id].append(feature)
            else:
                features.append(feature)

    # Merge multi-interval features
    for feature_id, intervals in multi_interval_features.items():
        if len(intervals) > 1:
            # Sort by start position
            intervals.sort(key=lambda x: x['start'])

            # For CDS features, exclude intervals that match intron coordinates
            if intervals[0]['type'] == 'CDS':
                valid_intervals = []
                for f in intervals:
                    coord = (f['start'], f['end'])
                    if coord not in intron_regions:
                        valid_intervals.append(f)
                    else:
                        print(f"  Excluding intron interval {coord} from CDS {intervals[0]['name']}")
                if valid_intervals:
                    intervals = valid_intervals

            # Create merged feature with multiple intervals
            merged = intervals[0].copy()
            merged['intervals'] = [(f['start'], f['end']) for f in intervals]
            features.append(merged)
        else:
            features.append(intervals[0])

    return features


def write_feature_table(features, seq_id, output_file):
    """Write features to NCBI feature table format."""
    with open(output_file, 'w') as f:
        # Write header
        f.write(f">Feature {seq_id}\n")

        # Sort features by start position
        features.sort(key=lambda x: x.get('intervals', [(x['start'], x['end'])])[0][0])

        for feature in features:
            feature_type = feature['type']
            strand = feature['strand']
            name = feature['name']

            # Map feature name if needed
            if feature_type in FEATURE_NAME_MAP:
                feature_type = FEATURE_NAME_MAP[feature_type]

            # Get intervals (single or multiple)
            if 'intervals' in feature:
                intervals = feature['intervals']
            else:
                intervals = [(feature['start'], feature['end'])]

            # For minus strand, reverse the order of start/end
            if strand == '-':
                intervals = [(end, start) for start, end in intervals]
                # Reverse interval order for - strand spliced features
                intervals = intervals[::-1]

            # Write first interval with feature type
            start, end = intervals[0]
            f.write(f"{start}\t{end}\t{feature_type}\n")

            # Write additional intervals (no feature type)
            for start, end in intervals[1:]:
                f.write(f"{start}\t{end}\n")

            # Write qualifiers
            gene_name = extract_gene_name(name)

            if feature_type == 'gene':
                if gene_name:
                    f.write(f"\t\t\tgene\t{gene_name}\n")

            elif feature_type == 'CDS':
                # CDS requires product qualifier
                if gene_name:
                    f.write(f"\t\t\tproduct\t{gene_name}p\n")
                    f.write(f"\t\t\tgene\t{gene_name}\n")
                else:
                    f.write(f"\t\t\tproduct\thypothetical protein\n")

                # Add codon_start if not starting at phase 0
                phase = feature.get('phase', '0')
                if phase and phase != '.' and phase != '0':
                    codon_start = int(phase) + 1
                    f.write(f"\t\t\tcodon_start\t{codon_start}\n")

            elif feature_type == 'mRNA':
                if gene_name:
                    f.write(f"\t\t\tproduct\t{gene_name}\n")
                    f.write(f"\t\t\tgene\t{gene_name}\n")

            elif feature_type == 'ncRNA':
                # ncRNA requires ncRNA_class qualifier
                if name:
                    f.write(f"\t\t\tproduct\t{name}\n")
                # Determine ncRNA class from name
                if 'SNR' in name.upper() or 'snR' in name:
                    f.write(f"\t\t\tncRNA_class\tsnoRNA\n")
                else:
                    f.write(f"\t\t\tncRNA_class\tother\n")
                if gene_name and gene_name != name:
                    f.write(f"\t\t\tgene\t{gene_name}\n")

            elif feature_type in ('tRNA', 'rRNA'):
                if name:
                    f.write(f"\t\t\tproduct\t{name}\n")
                if gene_name and gene_name != name:
                    f.write(f"\t\t\tgene\t{gene_name}\n")

            elif feature_type == 'regulatory':
                # regulatory requires regulatory_class qualifier
                # Try to extract class from name
                reg_class = 'other'
                name_lower = name.lower()
                if 'terminator' in name_lower:
                    reg_class = 'terminator'
                elif 'promoter' in name_lower:
                    reg_class = 'promoter'
                elif 'enhancer' in name_lower:
                    reg_class = 'enhancer'
                f.write(f"\t\t\tregulatory_class\t{reg_class}\n")
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'misc_binding':
                # misc_binding requires bound_moiety qualifier
                if name:
                    f.write(f"\t\t\tbound_moiety\t{name}\n")

            elif feature_type == 'modified_base':
                # modified_base requires mod_base qualifier
                f.write(f"\t\t\tmod_base\tother\n")
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'rep_origin':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'repeat_region':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'misc_feature':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'promoter':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'terminator':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")

            elif feature_type == 'primer_bind':
                if name:
                    f.write(f"\t\t\tnote\t{name}\n")


def get_sequence_id(fasta_file):
    """Extract sequence ID from FASTA header."""
    with open(fasta_file, 'r') as f:
        header = f.readline().strip()
        if header.startswith('>'):
            # Get first word after >
            seq_id = header[1:].split()[0]
            return seq_id
    return None


def get_sequence_length(fasta_file):
    """Get the sequence length from FASTA file."""
    seq_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq_length += len(line)
    return seq_length


def main():
    gff_file = 'mini-ess_XIV_circular.gff'
    fasta_file = 'mini-ess_XIV_circular.fasta'
    output_file = 'mini-ess_XIV_circular.tbl'

    # Get sequence ID from FASTA
    seq_id = get_sequence_id(fasta_file)
    if not seq_id:
        print("Error: Could not extract sequence ID from FASTA file")
        sys.exit(1)

    print(f"Sequence ID: {seq_id}")

    # Get sequence length
    seq_length = get_sequence_length(fasta_file)
    print(f"Sequence length: {seq_length} bp")

    # Parse GFF and write feature table
    features = parse_gff(gff_file)
    print(f"Parsed {len(features)} features")

    # Filter out features that extend beyond sequence boundaries or are too short
    valid_features = []
    skipped_coords = 0
    skipped_short = 0
    for f in features:
        if 'intervals' in f:
            max_coord = max(max(start, end) for start, end in f['intervals'])
            total_length = sum(abs(end - start) + 1 for start, end in f['intervals'])
        else:
            max_coord = max(f['start'], f['end'])
            total_length = abs(f['end'] - f['start']) + 1

        # Skip features beyond sequence boundaries
        if max_coord > seq_length:
            skipped_coords += 1
            print(f"  Skipping {f['type']} {f['name']}: coordinates exceed sequence length ({max_coord} > {seq_length})")
            continue

        # Skip CDSs that are too short (less than 3bp for a codon)
        if f['type'] == 'CDS' and total_length < 3:
            skipped_short += 1
            print(f"  Skipping {f['type']} {f['name']}: too short ({total_length}bp)")
            continue

        # Skip mRNA/gene that are too short (same as CDS issue)
        if f['type'] in ('mRNA', 'gene') and total_length < 3:
            skipped_short += 1
            print(f"  Skipping {f['type']} {f['name']}: too short ({total_length}bp)")
            continue

        valid_features.append(f)

    if skipped_coords:
        print(f"Skipped {skipped_coords} features with invalid coordinates")
    if skipped_short:
        print(f"Skipped {skipped_short} features that are too short")
    features = valid_features

    write_feature_table(features, seq_id, output_file)
    print(f"Feature table written to: {output_file}")

    # Print feature type counts
    type_counts = defaultdict(int)
    for f in features:
        type_counts[f['type']] += 1

    print("\nFeature counts:")
    for ftype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {ftype}: {count}")


if __name__ == '__main__':
    main()
