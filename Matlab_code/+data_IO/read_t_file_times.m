function spike_times = read_t_file_times(session_dir, t_file_names)
%READ_T_FILE_TIMES Read sorted Neuralynx/MClust T-file timestamps in seconds.

session_dir = char(string(session_dir));
if ischar(t_file_names) || isstring(t_file_names)
    t_file_names = cellstr(string(t_file_names));
end

spike_times = cell(numel(t_file_names), 1);
for file_index = 1:numel(t_file_names)
    t_file_name = char(string(t_file_names{file_index}));
    file_path = fullfile(session_dir, t_file_name);
    if exist(file_path, 'file') ~= 2
        error('data_IO:MissingTFile', 'T file not found: %s', file_path);
    end
    spike_times{file_index} = read_single_t_file(file_path);
end
end

function timestamps = read_single_t_file(file_path)

fid = fopen(file_path, 'r');
if fid < 0
    error('data_IO:UnreadableTFile', 'Could not open T file: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));

file_bytes = fread(fid, inf, '*uint8')';
header_marker = '%%ENDHEADER';
header_start = strfind(char(file_bytes), header_marker);
if isempty(header_start)
    error('data_IO:InvalidTFile', 'T-file header terminator was not found: %s', file_path);
end

payload_start = header_start(1) + numel(header_marker);
while payload_start <= numel(file_bytes) && ...
        (file_bytes(payload_start) == 10 || file_bytes(payload_start) == 13)
    payload_start = payload_start + 1;
end
payload = file_bytes(payload_start:end);
if isempty(payload)
    timestamps = zeros(0, 1);
    return
end

if mod(numel(payload), 8) == 0 && parse_unsigned(payload(1:8)) > double(intmax('uint32'))
    raw_timestamps = parse_timestamp_values(payload, 8);
    timestamps = raw_timestamps / 1e8;
elseif mod(numel(payload), 4) == 0
    raw_timestamps = parse_timestamp_values(payload, 4);
    timestamps = raw_timestamps / 1e4;
else
    error('data_IO:InvalidTFilePayload', 'T-file payload has an unsupported length: %s', file_path);
end
end

function values = parse_timestamp_values(payload, bytes_per_timestamp)

value_count = numel(payload) / bytes_per_timestamp;
values = zeros(value_count, 1);
for byte_index = 1:bytes_per_timestamp
    values = values * 256 + double(payload(byte_index:bytes_per_timestamp:end))';
end
end

function value = parse_unsigned(bytes)

value = 0;
for byte_index = 1:numel(bytes)
    value = value * 256 + double(bytes(byte_index));
end
end
