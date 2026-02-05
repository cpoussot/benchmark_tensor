function name = default_name(METHOD)

switch METHOD
    case 'mlf1' 
        name = 'A/G/P-V 2025 (A1)';
    case 'mlf2'
        name = 'A/G/P-V 2025 (A2)';
    case 'mdspack'
        name = 'MDSPACK v1.1.0';
    case 'kan1'
        name = 'P/P 2025';
    case 'paaa'
        name = 'C-R/B/G 2023';
    case 'paaalr'
        name = 'B/G 2025';
    case 'tensorflow'
        name = 'TensorFlow';
end
