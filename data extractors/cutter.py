file = 'new_cutsb.txt'
output = 'train_start2.exa'
before = 50
after = 50

with open(file) as file_handle, open(output, 'w') as output_handle:
    for line in file_handle:
        tokens = line.split()
        for it, token in enumerate(tokens[:1]):
            for i, char in enumerate(token):
                if char == 'P':
                    if i + after < len(token):
                        output_handle.write(token[i - before: i + after + 1].replace('P', '').lower() + '\n')
#                    else:
#                        to_take = i + after - len(token) + 1
#                        append = tokens[it + 1][:to_take]
#                        output_handle.write(token[i - before: i + after + 1].replace('P', '').lower() + append.lower() + '\n')
