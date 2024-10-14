def Reading(s):
    with open(s, 'r') as file:
        nums = []
        for line in file:
            num = line.split('<==')[0].strip()
            try:
                num = float(num)
                nums.append(num)
            except ValueError:
                print('ERROR ЧТЕНИЯ ФАЙЛА')
    return nums
