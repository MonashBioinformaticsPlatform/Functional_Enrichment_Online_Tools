document.addEventListener("DOMContentLoaded", function () {
    const headings = document.querySelectorAll('h1');
    headings.forEach(heading => {
        if (heading.textContent.includes("Day 1") || heading.textContent.includes("Day 2")) {
            heading.style.fontSize = "2.5em";
            heading.style.fontWeight = "bold";
            heading.style.textAlign = "center";
        }
    });
});
