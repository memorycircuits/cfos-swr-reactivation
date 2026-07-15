function [timestamps, samples, sampling_frequency] = read_ncs_downsampled(file_path, downsampling_factor)
%READ_NCS_DOWNSAMPLED Read samples and timestamps from a Neuralynx CSC file.
%   Timestamps are returned in seconds. Samples are returned at the requested
%   downsampled rate without loading the full CSC recording into memory.

if nargin < 2 || isempty(downsampling_factor)
    downsampling_factor = 1;
end
if ~isscalar(downsampling_factor) || downsampling_factor < 1 || ...
        downsampling_factor ~= floor(downsampling_factor)
    error('data_IO:InvalidDownsamplingFactor', ...
        'downsampling_factor must be a positive integer.');
end

file_path = char(string(file_path));
if exist(file_path, 'file') ~= 2
    error('data_IO:MissingCSCFile', 'CSC file not found: %s', file_path);
end

header_size = 16 * 1024;
record_size = 8 + 4 + 4 + 4 + 2 * 512;
file_info = dir(file_path);
data_size = file_info.bytes - header_size;
if data_size <= 0 || mod(data_size, record_size) ~= 0
    error('data_IO:InvalidCSCFile', 'Invalid Neuralynx CSC file structure: %s', file_path);
end

record_count = data_size / record_size;
samples_per_record = ceil(512 / downsampling_factor);
timestamps = zeros(record_count * samples_per_record, 1);
samples = zeros(record_count * samples_per_record, 1);
sampling_frequency = NaN;

fid = fopen(file_path, 'r', 'ieee-le');
if fid < 0
    error('data_IO:UnreadableCSCFile', 'Could not open CSC file: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));

fseek(fid, header_size, 'bof');
write_index = 0;
for record_index = 1:record_count
    record_timestamp = fread(fid, 1, 'uint64=>double');
    fread(fid, 1, 'uint32=>double');
    record_frequency = fread(fid, 1, 'uint32=>double');
    valid_sample_count = fread(fid, 1, 'uint32=>double');
    record_samples = fread(fid, 512, 'int16=>double');

    if isempty(record_timestamp) || isempty(record_frequency) || ...
            isempty(valid_sample_count) || numel(record_samples) ~= 512
        error('data_IO:TruncatedCSCFile', 'CSC record %d could not be read from %s.', record_index, file_path);
    end
    if record_frequency <= 0 || valid_sample_count > 512
        error('data_IO:InvalidCSCRecord', 'CSC record %d is invalid in %s.', record_index, file_path);
    end
    if isnan(sampling_frequency)
        sampling_frequency = record_frequency;
    elseif sampling_frequency ~= record_frequency
        error('data_IO:VariableSamplingFrequency', ...
            'CSC sampling frequency changes within %s.', file_path);
    end

    sample_indices = 1:downsampling_factor:valid_sample_count;
    sample_count = numel(sample_indices);
    destination_indices = write_index + (1:sample_count);
    timestamps(destination_indices) = record_timestamp / 1e6 + ...
        (sample_indices(:) - 1) / sampling_frequency;
    samples(destination_indices) = record_samples(sample_indices);
    write_index = write_index + sample_count;
end

timestamps = timestamps(1:write_index);
samples = samples(1:write_index);
sampling_frequency = sampling_frequency / downsampling_factor;
end
