## def header_parse
## parses Himawari HSD file header
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications:

def header_parse(header_bin, parse_blocks = [1,2,3,4,5,6,7,8,9,10,11]):
    import acolite as ac

    ## parse header
    ## get total number of header blocks (this should be fixed to 11)
    ret = ac.himawari.block_unpack(header_bin[3:5], [['total_number_of_header_blocks', 2, 'H']], strip = True)
    header_blocks = ret['total_number_of_header_blocks']

    ## get header block lengths (variable block lengths for certain headers)
    block_positions = []
    initial_position = 0
    for bi_ in range(header_blocks):
        bi = 1 + bi_
        ret = ac.himawari.block_unpack(header_bin[initial_position+1:initial_position+3],
                                       [['block_length', 2, 'H']], strip = True)
        block_length = ret['block_length']
        block_positions.append([initial_position, block_length])
        initial_position += block_length

    ## parse header data
    header_data = {}
    for bi in range(header_blocks):
        block = 1 + bi
        if block not in parse_blocks: continue

        ## subset header bin data
        bp = block_positions[bi]
        header_block = header_bin[bp[0]:bp[0]+bp[1]]

        ## block 1 - basic information
        if block == 1:
            b_bytes = [1, 2, 2, 1, 16, 16, 4, 2, 2, 8, 8, 8, 4, 4, 1, 1, 1, 1, 32, 128, 40]
            b_types = ['B', 'H', 'H', 'B', '%ds', '%ds', '%ds', '%ds', 'H', 'd', 'd', 'd', 'I', 'I', 'B', 'B', 'B', 'B', '%ds', '%ds', '%ds']
            b_names = ['header_block_number', 'block_length', 'total_number_of_header_blocks', 'byte_order']
            b_names += ['satellite_name', 'processing_center_name', 'observation_area', 'other_observation_information']
            b_names += ['observation_timeline', 'observation_start_time', 'observation_end_time', 'file_creation_time']
            b_names += ['total_header_length', 'total_data_length', 'quality_flag_1', 'quality_flag_2']
            b_names += ['quality_flag_3', 'quality_flag_4', 'file_format_version', 'file_name', 'spare']

        ## block 2 - data information
        elif block == 2:
            b_bytes = [1, 2, 2, 2, 2, 1, 40]
            b_types = ['B', 'H', 'H', 'H', 'H', 'B', '%ds']
            b_names = ['header_block_number', 'block_length', 'number_of_bits_per_pixel', 'number_of_columns']
            b_names += ['number_of_lines', 'compression_flag_for_data_block', 'spare']

        ## block 3 - projection
        elif block == 3:
            b_bytes = [1, 2, 8, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 2, 2, 40]
            b_types = ['B', 'H', 'd', 'I', 'I', 'I', 'I', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'H', 'H','%ds']
            b_names = ['header_block_number', 'block_length', 'sub_lon']
            b_names += ['CFAC', 'LFAC', 'COFF', 'LOFF']
            b_names += ['distance_from_earth_center', 'earth_equatorial_radius', 'earth_polar_radius']
            b_names += ['v11', 'v12', 'v13', 'v14', 'resampling_size', 'resampling_types', 'spare']

        ## block 4 - navigation
        elif block == 4:
            b_bytes = [1, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 40]
            b_types = ['B', 'H', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd','%ds']
            b_names = ['header_block_number', 'block_length', 'navigation_information_time']
            b_names += ['SSP_longitude', 'SSP_latitude', 'distance_from_earth_center']
            b_names += ['nadir_longitude', 'nadir_latitude']
            b_names += ['sun_position_x', 'sun_position_y', 'sun_position_z']
            b_names += ['moon_position_x', 'moon_position_y', 'moon_position_z']
            b_names += ['spare']

        ## block 5 - calibration
        elif block == 5:
            ## note that VNIR and TIR follow different header
            ## first get the band number
            ret = ac.himawari.block_unpack(header_block[3:5], [['band_number', 2, 'H']], strip = True)
            band_number = ret['band_number']

            ## use VNIR or TIR header
            if band_number in [1,2,3,4,5,6]:
                calibration_type = 'VNIR'
                ## block 5 - calibration (VNIR)
                b_bytes = [1, 2, 2, 8, 2, 2, 2, 8, 8, 8, 8, 8, 8, 80]
                b_types = ['B', 'H', 'H', 'd', 'H', 'H', 'H', 'd', 'd', 'd', 'd', 'd', 'd', '%ds']
                b_names = ['header_block_number', 'block_length', 'band_number']
                b_names += ['central_wavelength', 'valid_number_of_bits']
                b_names += ['count_value_of_error_pixels', 'count_value_of_pixels_outside_scan_area']
                b_names += ['slope_for_count_radiance_conversion', 'intercept_for_count_radiance_conversion']
                b_names += ['coefficient_for_radiance_to_albedo', 'update_time_of_the_values_12_13']
                b_names += ['calibrated_slope_for_count_radiance_conversion', 'calibrated_intercept_for_count_radiance_conversion']
                b_names += ['spare']
            elif band_number in [7,8,9,10,11,12,13,14,15,16]:
                calibration_type = 'TIR'
                ## block 5 - calibration (TIR) - not tested yet
                b_bytes = [1, 2, 2, 8, 2, 2, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 40]
                b_types = ['B', 'H', 'H', 'd', 'H', 'H', 'H', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', '%ds']
                b_names = ['header_block_number', 'block_length', 'band_number']
                b_names += ['central_wavelength', 'valid_number_of_bits']
                b_names += ['count_value_of_error_pixels', 'count_value_of_pixels_outside_scan_area']
                b_names += ['slope_for_count_radiance_conversion', 'intercept_for_count_radiance_conversion']
                b_names += ['correction_coefficient_c0', 'correction_coefficient_c1', 'correction_coefficient_c2']
                b_names += ['correction_coefficient_C0', 'correction_coefficient_C1', 'correction_coefficient_C2']
                b_names += ['speed_of_light_c', 'Planck_constant_h','Boltzmann_constant_k']
                b_names += ['spare']
            else:
                continue

        ## block 6 - intercalibration information
        elif block == 6:
            b_bytes = [1, 2, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 128, 56]
            b_types = ['B', 'H', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'f', 'f', '%ds', '%ds']
            b_names = ['header_block_number', 'block_length']
            b_names += ['GSICS_calibration_intercept', 'GSICS_calibration_slope', 'GSICS_calibration_quadratic']
            b_names += ['radiance_bias', 'radiance_bias_uncertainty', 'radiance_for_standard_scene']
            b_names += ['GSICS_validity_start_time', 'GSICS_validity_end_time']
            b_names += ['GSICS_radiance_validity_upper_limit', 'GSICS_radiance_validity_lower_limit']
            b_names += ['GSICS_file_name']
            b_names += ['spare']

        ## block 7 - segment information block
        elif block == 7:
            b_bytes = [1, 2, 1, 1, 2, 40]
            b_types = ['B', 'H', 'B', 'B', 'H', '%ds']
            b_names = ['header_block_number', 'block_length']
            b_names += ['total_number_of_segments', 'segment_sequence_number', 'first_line_number_of_image_sequence']
            b_names += ['spare']

        ## block 8 - navigation correction information block
        elif block == 8:
            ## note that elements 7-9 are repeated according to element 6
            ## first get the number of repeats
            ret = ac.himawari.block_unpack(header_block[19:21], [['number_of_correction_information_data_for_column_and_line_direction', 2, 'H']], strip = True)
            repeats = ret['number_of_correction_information_data_for_column_and_line_direction']

            ## parse data with repeating elements 7-9
            b_bytes = [1, 2, 4, 4, 8, 2]
            b_types = ['B', 'H', 'f', 'f', 'd', 'H']
            b_names = ['header_block_number', 'block_length']
            b_names += ['center_column_of_rotation', 'center_line_of_rotation', 'amount_of_rotational_correction']
            b_names += ['number_of_correction_information_data_for_column_and_line_direction']
            for r in range(repeats):
                b_bytes += [2, 4, 4]
                b_types += ['H', 'f', 'f']
                b_names += ['line_number_after_rotation_{}'.format(r+1),
                            'shift_for_column_direction_{}'.format(r+1),
                            'shift_for_line_direction_{}'.format(r+1)]
            ## add spare
            b_bytes += [40]
            b_types +=  ['%ds']
            b_names += ['spare']

        ## block 9 - observation time information block
        elif block == 9:
            ## note that elements 4-5 are repeated according to element 3
            ## first get the number of repeats
            ret = ac.himawari.block_unpack(header_block[3:5], [['number_of_observation_times', 2, 'H']], strip = True)
            repeats = ret['number_of_observation_times']

            b_bytes = [1, 2, 2]
            b_types = ['B', 'H', 'H']
            b_names = ['header_block_number', 'block_length']
            b_names += ['number_of_observation_times']

            ## add repeats
            for r in range(repeats):
                b_bytes += [2, 8]
                b_types += ['H', 'd']
                b_names += ['line_number_{}'.format(r+1),
                            'observation_time_{}'.format(r+1)]
            ## add spare
            b_bytes += [40]
            b_types +=  ['%ds']
            b_names += ['spare']

        ## block 10 - error information block
        elif block == 10:
            ## note that elements 4-5 are repeated according to element 3
            ## first get the number of repeats
            ret = ac.himawari.block_unpack(header_block[3:5], [['number_of_error_information_data', 2, 'H']], strip = True)
            repeats = ret['number_of_error_information_data']

            b_bytes = [1, 4, 2]
            b_types = ['B', 'I', 'H']
            b_names = ['header_block_number', 'block_length']
            b_names += ['number_of_error_information_data']

            ## add repeats
            for r in range(repeats):
                b_bytes += [2, 2]
                b_types += ['H', 'H']
                b_names += ['line_number_{}'.format(r+1),
                            'number_of_error_pixels_per_line_{}'.format(r+1)]
            ## add spare
            b_bytes += [40]
            b_types +=  ['%ds']
            b_names += ['spare']

        ## block 11 - spare block
        elif block == 11:
            b_bytes = [1, 2, 256]
            b_types = ['B', 'H', '%ds']
            b_names = ['header_block_number', 'block_length', 'spare']

        else:
            print('Block {} not configured'.format(block))
            continue

        ## parse block data
        header_data[block] = {}
        cur_el = 0
        cur_pos = 0
        while cur_pos < len(header_block)-1:
            cur_len = b_bytes[cur_el]
            ret = ac.himawari.block_unpack(header_block[cur_pos:cur_pos+cur_len], [[b_names[cur_el], b_bytes[cur_el], b_types[cur_el]]], strip = True)
            for k in ret.keys(): header_data[block][k] = ret[k]
            cur_el += 1
            cur_pos += cur_len

    return(header_data)
