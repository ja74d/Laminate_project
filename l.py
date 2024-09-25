import subprocess
from PIL import Image, ImageDraw, ImageFont

# Run the git log command and capture its output
result = subprocess.run(['git', 'log', '--oneline'], capture_output=True, text=True)
git_log_output = result.stdout

# Create an image with PIL
font = ImageFont.load_default()
# Estimate the size of the image
max_width = max(len(line) for line in git_log_output.split('\n'))
max_height = len(git_log_output.split('\n'))
image = Image.new('RGB', (max_width * 6, max_height * 10), color='black')
draw = ImageDraw.Draw(image)

# Draw the text on the image
draw.text((10, 10), git_log_output, font=font, fill='white')

# Save the image
image.save('git_log_output.png')

