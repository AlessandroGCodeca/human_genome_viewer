const puppeteer = require('puppeteer');

(async () => {
    const browser = await puppeteer.launch();
    const page = await browser.newPage();

    page.on('console', msg => console.log('PAGE LOG:', msg.text()));
    page.on('pageerror', error => console.log('PAGE ERROR:', error.message));
    page.on('requestfailed', request => console.log('REQUEST FAILED:', request.url(), request.failure().errorText));

    // Load the local HTML file
    const path = require('path');
    const filePath = `file://${path.resolve('test_ideogram.html')}`;
    console.log("Loading", filePath);

    await page.goto(filePath, { waitUntil: 'networkidle2' });
    const html = await page.evaluate(() => {
        return document.getElementById('ideo-container').innerHTML;
    });
    require('fs').writeFileSync('ideogram_dump.html', html);
    console.log("Wrote SVG to ideogram_dump.html");

    await browser.close();
})();
