function read_qtl(pathcasename, w2k_band_nums, atom_ind, char_inds, ef)
% read_qtl(pathcasename, w2k_band_nums, atom_ind, char_inds, ef)
%
% Reads .qtl and .energy files into a series of ascii files named pathcasename.band_N
% where N is the band number (not the WIEN2k band number)
%
% Also write a duplicate file to pathcasename.band_N_atom_M since cannot
% tell from file alone what atom is being examined
%
% Emulates Tony's Delphi program 'BandCharPlotter_v4'
%
% Format for the files from Tony's program is tab separated with no
% headers with each column as follows,
% kx, ky, kz, Energy, Total D, Dz2, Dxy, Dx2y2, Dxz+Dyz
% To get this, look at the first line of .qtl and choose the indexes for
% char_inds that correspond to the above (typically [11 7 8 9 10])
%
% pathcasename:    The filepath stem inclusing casename i.e. 'BaFe2P2_1e5/BaFe2P2' 
% w2k_band_nums:   The WIEN2k band numbers that you want to extract (n.b.
%                  will be renumbered to 1..n). If left empty, bands
%                  crossing Fermi energy will be used.
% atom_ind:        The atom number (as defined in case.struct, remove
%                  negative sign)
% char_inds:       The orbital characters to include in the files. Indexes
%                  are as laid out int the .qtl file header for a
%                  particular atom. Leave as [] to include all orbital
%                  characters.
% ef:              Fermi energy. If left empty, will use the Fermi energy
%                  from .qtl file (not usually a good idea!)

qtl_fn = [pathcasename '.qtl'];
disp(['QTL input file: ' qtl_fn]);

energy_fn = [pathcasename '.energy'];
disp(['Energy input file: ' energy_fn]);

qtl_fh = fopen(qtl_fn, 'r');
if qtl_fh == -1
    disp('ERROR: Cannot open QTL file!');
    return
end

w2k_band_nums = sort(w2k_band_nums);

% If no band numbers specified, will have to pre-read the QTL file to get
% the bands crossing the Fermi energy
if isempty(w2k_band_nums)
    disp('No band specified: Searching for bands crossing the Fermi energy ...');
    for dummy = 1:3
        line = fgetl(qtl_fh);
    end
    if isempty(ef)
        ef = str2double(line(57:66));
        disp(sprintf('  Fermi energy from QTL file: %fRy', ef));
    else
        disp(sprintf('  Fermi energy specified as: %fRy', ef));
    end
    % Skip to first band
    while ~strcmp(line(1:6), ' BAND:')
        line = fgetl(qtl_fh);
    end
    % Reads in all energies (even repeated) for each band, if they cross
    % the Fermi energy then add the band index
    w2k_band_num = 1;
    energies = [];
    while ischar(line)
        line = fgetl(qtl_fh);
        if ~ischar(line) || strcmp(line(1:6), ' BAND:')
            if (max(energies) > ef) && (min(energies) < ef)
                disp(sprintf('  WIEN2k band %d crosses Fermi energy', w2k_band_num));
                w2k_band_nums = [w2k_band_nums w2k_band_num];
            end
            w2k_band_num = w2k_band_num + 1;
            energies = [];
        else
            energies = [energies str2double(line(1:9))];
        end
    end
    % Reopen the QTL file for use
    fclose(qtl_fh);
    qtl_fh = fopen(qtl_fn, 'r');
end


disp('Reading QTL file ...');
% Read in the number of atoms from the fourth line '... NAT=...'
for dummy = 1:4
    line = fgetl(qtl_fh);
end
num_atoms = str2num(line(36:39));
% Read character data from the appropriate JATOM line
for atom_line = 1:atom_ind
    line = fgetl(qtl_fh);
end
char_names = textscan(strtrim(line(32:end)), '%s', 'Delimiter', ',');
char_names = char_names{1};
disp('The following orbital characters were found:');
for ind = 1:length(char_names)
    disp(sprintf('%2d: %s', ind, char_names{ind}));
end
disp('Selected the following:');
if isempty(char_inds)
    disp('ALL');
else
    disp(char_inds);
end
% Now read in the character data for each band
disp('Reading QTL file ...');
qtl_data = {};
for w2k_band_num = w2k_band_nums
    % Fast forward to the first band in the list
    while ~strcmp(line(1:10), sprintf(' BAND:%4d', w2k_band_num));
        line = fgetl(qtl_fh);
    end
    line = fgetl(qtl_fh);
    % Read lines until next band reached
    line_num = 1;
    lines_per_kpt = num_atoms + 1;
    qtl_data{end + 1} = [];
    while ~strcmp(line(1:6), ' BAND:');
        % Read in data only from relavent atom
        if mod(line_num - atom_ind, lines_per_kpt) == 0
            % Remove the weird extra spaces in the columns
            line = [line(14:21) line(25:end)];
            % Rest of line is made up of width 8 floats, one for each
            % character
            qtl_data{end} = [qtl_data{end}; sscanf(line, '%8f')'];
        end
        line = fgetl(qtl_fh);
        % Exit loop if EOF
        if ~ischar(line)
            break;
        end
        line_num = line_num + 1;
    end
end
fclose(qtl_fh);
num_kpts = size(qtl_data{1}, 1);


disp('Reading ENERGY file ...');

energy_fh = fopen(energy_fn, 'r');
if energy_fh == -1
    disp('ERROR: Cannot open ENERGY file!');
    return
end

% Skip 6 lines of headers to get to first k-point line
for dummy = 1:7
    line = fgetl(energy_fh);
end
energy_data = [];
skip_lines = 0;
lines_til_kpt = 0;
while ischar(line)
    if lines_til_kpt == 0
        % Is a k-point line containing the kx, ky, kz values and the number
        % of lines before the next k-point line
        line_data = sscanf(line, '%19f', 3);
        line_data = [line_data' NaN(1, length(w2k_band_nums))];
        lines_til_kpt = str2double(line(74:79)) + 1;
    else
        % Is a band energy line. Only read in the data if asked for
        band_num = str2double(line(1:12));
        band_ind = find(band_num == w2k_band_nums);
        if ~isempty(band_ind)
            line_data(3 + band_ind) = str2double(line(13:end));
        end
    end
    line = fgetl(energy_fh);
    lines_til_kpt = lines_til_kpt - 1;
    if lines_til_kpt == 0
        energy_data = [energy_data; line_data];
    end
end
fclose(energy_fh);

% Write the data to files
for band_num = 1:length(w2k_band_nums)
    if isempty(char_inds)
        char_data = qtl_data{band_num};
    else
        char_data = qtl_data{band_num}(:, char_inds);
    end
    out_data = [energy_data(:,1:3) energy_data(:,3+band_num) char_data];
    out_fn = [pathcasename '.band_' num2str(band_num)];
    save(out_fn, 'out_data', '-ascii', '-tabs');
    disp(['Written to: ' out_fn]);
    out_fn = [pathcasename '.band' num2str(band_num)  '_atom' num2str(atom_ind)];
    save(out_fn, 'out_data', '-ascii', '-tabs');
    disp(['Also written to: ' out_fn]);
end

disp('read_qtl done.');

end



