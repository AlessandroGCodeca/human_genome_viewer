from google import genai
from google.genai import types

def get_ai_response(gene_id: str, query: str, api_key: str = None) -> str:
    """
    Use Google's Gemini to explain gene function.
    """
    if not api_key:
        return "⚠️ Please enter a Google Gemini API Key in the text box above to use the AI Assistant."
        
    try:
        client = genai.Client(api_key=api_key)
        
        system_instruction = (
            f"You are an expert bioinformatics AI assistant embedded in a Human Genome Viewer application. "
            f"The user is currently analyzing the gene or accession ID: {gene_id}. "
            f"Answer the user's question clearly, scientifically, and concisely. "
            f"If the question is completely unrelated to genetics, biology, or the application, gently steer them back."
        )
        
        response = client.models.generate_content(
            model='gemini-2.5-flash',
            contents=query,
            config=types.GenerateContentConfig(
                system_instruction=system_instruction,
                temperature=0.3,
            )
        )
        return response.text
        
    except Exception as e:
        return f"❌ An error occurred while communicating with Gemini API: {str(e)}"
