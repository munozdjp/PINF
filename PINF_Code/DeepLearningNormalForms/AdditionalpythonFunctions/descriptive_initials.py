# Contents of descriptive_initials.py

import matplotlib.pyplot as plt



def get_descriptive_initials(filename):
    descriptor_map = {
        'XSadd': 'S',
        'XPitch': 'P',
        'XTrans': 'T',
        'XHopf': 'H',
        'XHopfSaddBig': 'SH',
        'XSaddPitchfork': 'SP',
        'XSaddTranscritical': 'ST',
        'XPitchforkTranscritical': 'PT',
        'XPitchforkHopf': 'PH',
        'XTranscriticalHopf': 'TH'
    }

    terms = filename.replace('Data_And_labels_', '').replace('.mat', '').replace('_Final_', ':').split('_')

    # Expand the terms by splitting further on ':'
    expanded_terms = []
    for term in terms:
        expanded_terms.extend(term.split(':'))

    initials = [descriptor_map[term] for term in expanded_terms if term in descriptor_map]
    short_title = '-'.join(initials)
    return short_title
